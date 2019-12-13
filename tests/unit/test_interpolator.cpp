#include <gtest/gtest.h>
#include <interpolation.h>
#include <gsl_interpolation_interface.h>

// ========== Test interpolator on linear grid ==============================================================================

class lin_interpolatorTest : public ::testing::Test,
			     public scam::utl::interpolator< scam::utl::gsl_lin_interp > {  
  
protected:

  lin_interpolatorTest () :
    interpolator< scam::utl::gsl_lin_interp >{ scam::utl::lin_vector( 20, 0.5, 10. ), scam::utl::lin_vector( 20, 0.5, 10. ) }
  {

  };

  virtual void SetUp() {

    my_class_test_ = this;

  }

  lin_interpolatorTest * my_class_test_;

};

TEST_F( lin_interpolatorTest, Constructor ) {
  EXPECT_DOUBLE_EQ( 0.5,
		    ( *my_class_test_ )( 0.5 ) );
  EXPECT_NEAR( 1.75,
	       ( *my_class_test_ )( 1.75 ),
	       1.e-4 );
  EXPECT_NEAR( 1.5,
	       my_class_test_->integrate( 1., 2. ),
	       1.e-4 );
}

// ========== Test interpolator on logarithmic grid =========================================================================

class log_interpolatorTest : public ::testing::Test,
			     public scam::utl::interpolator< scam::utl::gsl_log_interp > {  
  
protected:

  log_interpolatorTest () :
    interpolator< scam::utl::gsl_log_interp >{ scam::utl::log_vector( 20, 0.5, 10. ), scam::utl::log_vector( 20, 0.5, 10. ) }
  {

  };

  virtual void SetUp() {

    my_class_test_ = this;

  }

  log_interpolatorTest * my_class_test_;

};

TEST_F( log_interpolatorTest, Constructor ) {
  EXPECT_DOUBLE_EQ( 0.5,
		    ( *my_class_test_ )( 0.5 ) );
  EXPECT_NEAR( 1.75,
	       ( *my_class_test_ )( 1.75 ),
	       1.e-4 );
  EXPECT_NEAR( 1.5,
	       my_class_test_->integrate( 1., 2. ),
	       1.e-4 );
}

// ==========================================================================================================================
