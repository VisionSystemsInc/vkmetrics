#include <testlib/testlib_register.h>

DECLARE( test_binary_regions );
DECLARE( test_containment_tree );
DECLARE( test_obj_io );
DECLARE( test_ground_truth );
DECLARE( test_metrics );

void
register_tests()
{
  REGISTER( test_binary_regions );
  REGISTER( test_containment_tree );
  REGISTER( test_obj_io );
  REGISTER( test_ground_truth );
  REGISTER( test_metrics );
}

DEFINE_MAIN;
