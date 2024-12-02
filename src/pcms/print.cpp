#include "print.h"

namespace pcms {

  FILE* pcms_stdout = stdout;
  FILE* pcms_stderr = stderr;

  FILE* getStdout() { return pcms_stdout; }
  FILE* getStderr() { return pcms_stderr; }

  void setStdout(FILE* out) {
    assert(out != NULL);
    pcms_stdout = out;
  }

  void setStderr(FILE* err) {
    assert(err != NULL);
    pcms_stderr = err;
  }
}