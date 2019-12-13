#include <gtest/gtest.h>
#include <cosmology.h>
// #include <gsl_interpolation_interface.h>

// ========== Test interpolator on linear grid ==============================================================================

class cosmologyTest : public ::testing::Test,
		      public scam::cosmology {  
  
protected:

  cosmologyTest () :
    cosmology {}
  {

  };

  virtual void SetUp() {

    my_class_test_ = this;

  }

  cosmologyTest * my_class_test_;

};

TEST_F( cosmologyTest, Constructor ) {
  EXPECT_DOUBLE_EQ( 0.,
		    my_class_test_->d_C( 1.e-7 ) );
  // EXPECT_NEAR( 1.75,
  // 	       ( *my_class_test_ )( 1.75 ),
  // 	       1.e-4 );
  // EXPECT_NEAR( 1.5,
  // 	       my_class_test_->integrate( 1., 2. ),
  // 	       1.e-4 );
}
