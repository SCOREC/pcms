#include <catch2/catch_test_macros.hpp>
#include "pcms/utility/assert.h"
#include "pcms/utility/print.h"
#include <iostream>

int raise_error(int code)
{
  if (code) {
    throw pcms::pcms_error("Test exception - Raising error for testing");
    return 1;
  }
  return 0;
}

TEST_CASE("pcms error handling test with try-catch")
{
  bool exception_caught = false;
  bool correct_exception_type = false;

  try {
    raise_error(1);
    FAIL("Expected exception was not thrown");
  } catch (const pcms::pcms_error& e) {
    exception_caught = true;
    correct_exception_type = true;
    std::cout << "Caught pcms_error: " << e.what() << std::endl;
  }

  REQUIRE(exception_caught);
  REQUIRE(correct_exception_type);
}
