#include <iostream>

#include <helib/helib.h>
#include <helib/binaryArith.h>
#include <helib/intraSlot.h>
#include <helib/binaryCompare.h>

int main(int argc, char* argv[])
{
  long p = 2;
  // Cyclotomic polynomial - defines phi(m).
  long m = 4095;
  // Hensel lifting (default = 1).
  long r = 1;
  // Number of bits of the modulus chain.
  long bits = 500;
  // Number of columns of Key-Switching matrix (typically 2 or 3).
  long c = 2;
  // Factorisation of m required for bootstrapping.
  std::vector<long> mvec = {7, 5, 9, 13};
  // Generating set of Zm* group.
  std::vector<long> gens = {2341, 3277, 911};
  // Orders of the previous generators.
  std::vector<long> ords = {6, 4, 6};

  std::cout << "Initialising context object..." << std::endl;
  // Initialize the context.
  helib::Context context(m, p, r, gens, ords);

  // Modify the context, adding primes to the modulus chain.
  std::cout << "Building modulus chain..." << std::endl;
  buildModChain(context, bits, c);

  // Make bootstrappable.
  context.makeBootstrappable(
      helib::convert<NTL::Vec<long>, std::vector<long>>(mvec));

  // Print the context.
  context.zMStar.printout();
  std::cout << std::endl;

  // Print the security level.
  std::cout << "Security: " << context.securityLevel() << std::endl;

  // Secret key management.
  std::cout << "Creating secret key..." << std::endl;
  // Create a secret key associated with the context.
  helib::SecKey secret_key(context);
  // Generate the secret key.
  secret_key.GenSecKey();

  // Generate bootstrapping data.
  secret_key.genRecryptData();

  // Public key management.
  // Set the secret key (upcast: SecKey is a subclass of PubKey).
  const helib::PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context.
  const helib::EncryptedArray& ea = *(context.ea);

  // Build the unpack slot encoding.
  std::vector<helib::zzX> unpackSlotEncoding;
  buildUnpackSlotEncoding(unpackSlotEncoding, ea);

  // Get the number of slot (phi(m)).
  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;


  long bitSize = 16;
  long outSize = 2 * bitSize;

  long P_x = 20; 
  long P_y = 15;
  long R_x = 10; 
  long R_y = 10; 

  std::cout << "Pre-encryption data:" << std::endl;
  std::cout << "P_x = " << P_x << "  P_y = " << P_y << std::endl;
  std::cout << "R_x = " << R_x << "  R_y = " << R_y << std::endl;

  // Use a scratch ciphertext to populate vectors.
  helib::Ctxt scratch(public_key);
  std::vector<helib::Ctxt> encrypted_P_x(bitSize, scratch);
  std::vector<helib::Ctxt> encrypted_P_y(bitSize, scratch);
  std::vector<helib::Ctxt> encrypted_R_x(bitSize, scratch);
  std::vector<helib::Ctxt> encrypted_R_y(bitSize, scratch);
  // Encrypt the data in binary representation.
  for (long i = 0; i < bitSize; ++i) {
    std::vector<long> P_x_vec(ea.size());
    std::vector<long> P_y_vec(ea.size());
    std::vector<long> R_x_vec(ea.size());
    std::vector<long> R_y_vec(ea.size());
    for (auto& slot : P_x_vec)
      slot = (P_x >> i) & 1;
    for (auto& slot : P_y_vec)
      slot = (P_y >> i) & 1;
    for (auto& slot : R_x_vec)
      slot = (R_x >> i) & 1;
    for (auto& slot : R_y_vec)
      slot = (R_y >> i) & 1;
    ea.encrypt(encrypted_P_x[i], public_key, P_x_vec);
    ea.encrypt(encrypted_P_y[i], public_key, P_y_vec);
    ea.encrypt(encrypted_R_x[i], public_key, R_x_vec);
    ea.encrypt(encrypted_R_y[i], public_key, R_y_vec);
  }

// ************************************************************************
  // Compute delta_x = E(P_x) - E(R_x)
  std::vector<helib::Ctxt> encrypted_delta_x(bitSize, scratch);
  helib::CtPtrs_vectorCt wrapper_delta_x(encrypted_delta_x);
  helib::subtractBinary(
      wrapper_delta_x,
      helib::CtPtrs_vectorCt(encrypted_P_x),
      helib::CtPtrs_vectorCt(encrypted_R_x),
      &unpackSlotEncoding);

  // Compute delta_y = E(P_y) - E(R_y)
  std::vector<helib::Ctxt> encrypted_delta_y(bitSize, scratch);
  helib::CtPtrs_vectorCt wrapper_delta_y(encrypted_delta_y);
  helib::subtractBinary(
      wrapper_delta_y,
      helib::CtPtrs_vectorCt(encrypted_P_y),
      helib::CtPtrs_vectorCt(encrypted_R_y),
      &unpackSlotEncoding);

  // Compute the square of delta_x and delta_y
  std::vector<helib::Ctxt> encrypted_square_delta_x;
  helib::CtPtrs_vectorCt wrapper_square_delta_x(encrypted_square_delta_x);
  helib::multTwoNumbers(
      wrapper_square_delta_x,
      wrapper_delta_x,
      wrapper_delta_x,
      false, 
      outSize, 
      &unpackSlotEncoding);
  std::vector<helib::Ctxt> encrypted_square_delta_y;
  helib::CtPtrs_vectorCt wrapper_square_delta_y(encrypted_square_delta_y);
  helib::multTwoNumbers(
      wrapper_square_delta_y,
      wrapper_delta_y,
      wrapper_delta_y,
      false, 
      outSize, 
      &unpackSlotEncoding);

  // Sum of two squared numbers
  std::vector<helib::Ctxt> encrypted_result;
  helib::CtPtrs_vectorCt result_wrapper(encrypted_result);
  helib::addTwoNumbers(
      result_wrapper,
      wrapper_square_delta_x,
      wrapper_square_delta_y,
      false,
      &unpackSlotEncoding);

  // Decrypt and print the result.
  std::vector<long> decrypted_result;
  helib::decryptBinaryNums(decrypted_result, result_wrapper, secret_key, ea);
  std::cout << "D(d(E(P),E(R)) = " << decrypted_result.back() << std::endl;

  std::cout << "d(P,R) = " << (P_x-R_x)*(P_x-R_x)+(P_y-R_y)*(P_y-R_y) << std::endl;

  return 0;
}
