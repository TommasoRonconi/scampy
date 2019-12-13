#include <gtest/gtest.h>
#include <harikane16_p.h>
#include <tinker10_p.h>

class harikane16_pTest: public ::testing::Test, public scam::harikane16_p {
  
protected:

  harikane16_pTest () :
    harikane16_p{ 0.5, 1.e+11, 1., 1.e+11, 1.e+11, 1. }
  {
    
  };
  
  virtual void SetUp() {
    
    my_class_test_ = this;
    
  }
  
  harikane16_pTest * my_class_test_;
  
};

TEST_F( harikane16_pTest, Constructor ) {
  EXPECT_DOUBLE_EQ( 0.9213503964748575,
		    my_class_test_->Ncen( 1.e+12 ) );
  EXPECT_DOUBLE_EQ( 8.2921535682737169,
		    my_class_test_->Nsat( 1.e+12 ) );
}

