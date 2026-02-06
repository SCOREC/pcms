#include <catch2/catch_test_macros.hpp>
#include <pcms/print.h>
#include <Kokkos_Core.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <filesystem>
#include <mpi.h>

// Helper function to create unique temp file path
std::string create_temp_filename() {
  return (std::filesystem::temp_directory_path() /
         ("pcms_test_" + std::to_string(std::rand()) + ".txt")).string();
}

TEST_CASE("printInfo on host", "[print]") {
  // Create temporary file for output
  std::string temp_path = create_temp_filename();

  // Need FILE* for PCMS API, but will use C++ streams to read
  FILE* temp_file = std::fopen(temp_path.c_str(), "w+");
  REQUIRE(temp_file != nullptr);

  // Redirect stdout to our test file
  FILE* original_stdout = pcms::getStdout();
  pcms::setStdout(temp_file);

  SECTION("printInfo writes to redirected stdout") {
    const char* test_message = "Test message from rank %d\n";
    int rank = 42;

    // Print to redirected stdout
    pcms::printInfo(test_message, rank);

    // Flush and close to ensure data is written
    std::fflush(temp_file);

    // Read back using C++ streams
    std::ifstream input(temp_path);
    REQUIRE(input.is_open());

    std::string line;
    std::getline(input, line);

    REQUIRE(line == "Test message from rank 42");
    input.close();
  }

  SECTION("multiple printInfo calls accumulate output") {
    pcms::printInfo("Line 1: %d\n", 1);
    pcms::printInfo("Line 2: %d\n", 2);
    pcms::printInfo("Line 3: %d\n", 3);

    std::fflush(temp_file);

    // Read back using C++ streams
    std::ifstream input(temp_path);
    REQUIRE(input.is_open());

    std::string line;
    int line_count = 0;
    while (std::getline(input, line)) {
      line_count++;
    }

    REQUIRE(line_count == 3);
    input.close();
  }

  // Restore original stdout and cleanup
  pcms::setStdout(original_stdout);
  std::fclose(temp_file);
  std::filesystem::remove(temp_path);
}

TEST_CASE("printInfo in device code", "[print][kokkos]") {
  SECTION("printInfo compiles and runs in Kokkos parallel kernel") {
    const int n = 10;
    int success_count = 0;

    // This test verifies that printInfo can be called from device code
    // without compilation errors. On device, it should be a no-op.
    Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(int i, int& local_count) {
        pcms::printInfo("Device iteration %d\n", i); //no-op
        local_count++;
      },
      success_count
    );

    // Verify the kernel executed for all iterations
    REQUIRE(success_count == n);
  }
}

TEST_CASE("printInfo works on host in GPU builds", "[print]") {
  // This test ensures that even when compiled for GPU execution,
  // printInfo still works on the host side

  std::string temp_path = create_temp_filename();
  FILE* temp_file = std::fopen(temp_path.c_str(), "w+");
  REQUIRE(temp_file != nullptr);

  FILE* original_stdout = pcms::getStdout();
  pcms::setStdout(temp_file);

  SECTION("host printInfo works before kernel launch") {
    pcms::printInfo("Before kernel: %d\n", 123);
    std::fflush(temp_file);

    std::ifstream input(temp_path);
    REQUIRE(input.is_open());

    std::string line;
    std::getline(input, line);

    REQUIRE(line == "Before kernel: 123");
    input.close();
  }

  SECTION("host printInfo works after kernel launch") {
    // Launch a kernel (which may call printInfo as no-op)
    Kokkos::parallel_for(10, KOKKOS_LAMBDA(int i) {
      pcms::printInfo("Device: %d\n", i);  // No-op on device
    });
    Kokkos::fence();

    // Host printInfo should still work
    pcms::printInfo("After kernel: %d\n", 456);
    std::fflush(temp_file);

    std::ifstream input(temp_path);
    REQUIRE(input.is_open());

    std::string line;
    std::getline(input, line);

    // Should only see the host message, not the 10 device messages
    REQUIRE(line == "After kernel: 456");

    // Verify there's only one line (no device output)
    int line_count = 1;
    while (std::getline(input, line)) {
      line_count++;
    }
    REQUIRE(line_count == 1);
    input.close();
  }

  pcms::setStdout(original_stdout);
  std::fclose(temp_file);
  std::filesystem::remove(temp_path);
}

TEST_CASE("printError redirects to stderr", "[print]") {
  std::string temp_path = create_temp_filename();
  FILE* temp_file = std::fopen(temp_path.c_str(), "w+");
  REQUIRE(temp_file != nullptr);

  FILE* original_stderr = pcms::getStderr();
  pcms::setStderr(temp_file);

  pcms::printError("Error: code %d\n", 404);
  std::fflush(temp_file);

  std::ifstream input(temp_path);
  REQUIRE(input.is_open());

  std::string line;
  std::getline(input, line);

  REQUIRE(line.find("Error: code 404") != std::string::npos);
  input.close();

  pcms::setStderr(original_stderr);
  std::fclose(temp_file);
  std::filesystem::remove(temp_path);
}

TEST_CASE("rank-specific file output example", "[print]") {
  // Example showing how to redirect each rank's output to a separate file
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream ss;
  ss << std::filesystem::temp_directory_path().string()
     << "/output_rank_" << rank << ".txt";
  std::string rank_file_path = ss.str();

  FILE* rank_file = std::fopen(rank_file_path.c_str(), "w");
  REQUIRE(rank_file != nullptr);

  FILE* original_stdout = pcms::getStdout();
  pcms::setStdout(rank_file);

  // All printInfo calls will now go to rank-specific file
  pcms::printInfo("Message from rank %d\n", rank);
  pcms::printInfo("Another message from rank %d\n", rank);

  std::fflush(rank_file);

  // Verify output
  std::ifstream input(rank_file_path);
  REQUIRE(input.is_open());

  std::stringstream buffer;
  buffer << input.rdbuf();
  std::string contents = buffer.str();

  REQUIRE(contents.find("Message from rank") != std::string::npos);
  input.close();

  // Cleanup
  pcms::setStdout(original_stdout);
  std::fclose(rank_file);
  std::filesystem::remove(rank_file_path);
}
