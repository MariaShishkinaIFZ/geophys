#include <iostream>
#include <stdlib.h>
#include <values.h>
#include "gsterror.h"

using namespace std;

void GSTerror(const char *s)
{
  cerr << s;
  exit(EXIT_FAILURE);
}

void GSTwarning(const char *s)
{
  cerr << s;
}
