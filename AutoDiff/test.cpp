#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "AutoDiff.h"
#include <complex>
#include "BlackScholesAutoDiff.h"

TEST_CASE("Test polynomial", "[Functional]"){
    AutoDiff<double> myAutoDiff(5.0, 1.0);
    auto mySquare=[&](const auto& myX){
        return myX*myX;
    };
    auto myResult=mySquare(myAutoDiff);


    REQUIRE(myResult.getStandard()==Approx(25.0));
    REQUIRE(myResult.getDual()==Approx(10.0));
}
TEST_CASE("Test <<", "[Functional]"){
    AutoDiff<double> myAutoDiff(5.0, 1.0);
    std::cout<<myAutoDiff<<std::endl;
    

}
TEST_CASE("Test exponential", "[Functional]"){
    AutoDiff<double> myAutoDiff(5.0, 1.0);
    auto myExp=[&](const auto& myX){
        return exp(myX);
    };
    auto myResult=myExp(myAutoDiff);


    REQUIRE(myResult.getStandard()==Approx(exp(5.0)));
    REQUIRE(myResult.getDual()==Approx(exp(5.0)));
}
TEST_CASE("Test abs", "[Functional]"){
    AutoDiff<double> myAutoDiff(-5.0, 1.0);
    auto myResult=abs(myAutoDiff);
    REQUIRE(myResult==Approx(5.0));

}
TEST_CASE("Test fabs", "[Functional]"){
    AutoDiff<double> myAutoDiff(-5.0, 1.0);
    auto myResult=fabs(myAutoDiff);
    REQUIRE(myResult==Approx(5.0));

}
TEST_CASE("Test logistic", "[Functional]"){
    AutoDiff<double> myAutoDiff(.5, 1.0);
    auto myF=[](const auto& val){
        return 1.0/(1.0+exp(-val));
    };
    auto myResult=myF(myAutoDiff);
    auto expected=myF(.5);
    REQUIRE(myResult.getStandard()==Approx(expected));
    REQUIRE(myResult.getDual()==Approx(expected*(1.0-expected)));
}

TEST_CASE("Test complex", "[Functional]"){
    AutoDiff<std::complex<double> > myAutoDiff(std::complex<double>(0.0, 4.0), std::complex<double>(0.0, 1.0));
    auto myExp=[&](const auto& myX){
        return exp(myX);
    };
    auto myResult=myExp(myAutoDiff);


    REQUIRE(myResult.getStandard().real()==Approx(cos(4.0)));
    REQUIRE(myResult.getStandard().imag()==Approx(sin(4.0)));
    REQUIRE(myResult.getDual().real()==Approx(-sin(4.0)));
    REQUIRE(myResult.getDual().imag()==Approx(cos(4.0)));
}

TEST_CASE("Test Delta", "[BS]"){
    const auto S0=AutoDiff<double>(50, 1);
    const auto r=.05;
    const auto t=1.0;
    const auto discount=exp(-r*t);
    const auto k=50;
    const auto sigma=.3*sqrt(t);

    double s=sqrt(2.0);
    auto d1=log(50.0/(discount*k))/(sigma)+sigma*.5;
    REQUIRE(autodiffonly::BSCall(S0, discount, k, sigma).getDual()==Approx(.5*erf(d1/s)+.5));
}

