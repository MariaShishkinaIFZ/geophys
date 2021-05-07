#include <math.h>

#include "gstdefs.h"
#include "gstmath.h"

bool GSTfloatsEqual(float x1, float x2)
{
  if (x1+x2>GST_FLOAT_ACCURACY)
    return (0.5*fabs(x1-x2)/(x1+x2) <= GST_FLOAT_ACCURACY);
  else return true;
}

bool GSTdoublesEqual(double x1, double x2)
{
  if (x1+x2>GST_DOUBLE_ACCURACY)
    return (0.5*fabs(x1-x2)/(x1+x2) <= GST_DOUBLE_ACCURACY);
  else return true;
}

bool GSTnumsEqual(float x1, float x2)
{
  if (x1+x2>GST_FLOAT_ACCURACY)
    return (0.5*fabs(x1-x2)/(x1+x2) <= GST_FLOAT_ACCURACY);
  else return true;
}

bool GSTnumsEqual(double x1, double x2)
{
  if (x1+x2>GST_DOUBLE_ACCURACY)
    return (0.5*fabs(x1-x2)/(x1+x2) <= GST_DOUBLE_ACCURACY);
  else return true;
}

bool GSTnumsEqual(int x1, int x2)
{
  if (x1!=x2) return false;
  else return true;
}

bool GSTnonZero(float x1)
{
  if (x1>GST_FLOAT_ACCURACY) return true;
  else return false;
}

bool GSTnonZero(double x1)
{
  if (x1>GST_DOUBLE_ACCURACY) return true;
  else return false;
}

bool GSTnonZero(int x1)
{
  if (x1!=0) return true;
  else return false;
}
