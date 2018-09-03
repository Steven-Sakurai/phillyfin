#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef AUTODIFF_H
#define AUTODIFF_H


#include <cmath>
template<class T>
class AutoDiff {
    private:
        T standard;
        T dual;

    public:
        AutoDiff(const T& stnd, const T& dl){
            standard=stnd;
            dual=dl;
        }
        AutoDiff(){
            standard=0;
            dual=0;
        }

        AutoDiff(const T& stnd){
            standard=stnd;
            dual=0;
        }
        T getStandard() const{
            return standard;
        }
        T getDual() const{
            return dual;
        }
        AutoDiff<T> recipricol() const{
            return AutoDiff<T>(1.0/standard, -dual/(standard*standard));
        }

        void setDual(const T& val){
            dual=val;
        }
        void setStandard(const T& val){
            standard=val;
        }
        AutoDiff<T>& operator=(const AutoDiff<T>& right){
            standard=right.getStandard();
            dual=right.getDual();
            return *this;
        }
        AutoDiff<T>& operator=(const T& x){
            standard=x;
            dual=0.0;
            return *this;
        }
        template<typename W>
        auto add(const AutoDiff<W> &val) const{
            return AutoDiff(standard+val.getStandard(), dual+val.getDual());
        }
        template< typename W>
        auto add(const W& val) const{
            return AutoDiff(standard+val, dual);
        }
        template<typename W>
        auto subtract(const AutoDiff<W> &val) const{
            return AutoDiff(standard-val.getStandard(), dual-val.getDual());
        }
        template<typename W>
        auto subtract(const W& val) const{
            return AutoDiff(standard-val, dual);
        }
        template<typename W>
        auto multiply(const AutoDiff<W> &val) const{
            return AutoDiff(standard*val.getStandard(), standard*val.getDual()+dual*val.getStandard());
        }
        template<typename W>
        auto multiply(const W& val) const{
            return AutoDiff(standard*val, dual*val);
        }
        template<typename W>
        void operator+=(const AutoDiff<W> &right){
            standard+=right.getStandard();
            dual+=right.getDual();
        }
        friend std::ostream &operator<<(std::ostream& out, const AutoDiff<T>& objToPrint){
            out<<"standard: "<<objToPrint.getStandard()<<", dual: "<<objToPrint.getDual();
            return out;
        }
        
};





template< typename W>
auto sin(const AutoDiff<W> &val){
    return AutoDiff<W>(sin(val.getStandard()), cos(val.getStandard())*val.getDual());
}
template<typename W>
auto cos(const AutoDiff<W> &val){
    return AutoDiff<W>(cos(val.getStandard()), -sin(val.getStandard())*val.getDual());
}
template<typename W>
auto exp(const AutoDiff<W> &val){
    return AutoDiff<W>(exp(val.getStandard()), exp(val.getStandard())*val.getDual());
}
template<typename W>
auto log(const AutoDiff<W> &val){
    return AutoDiff<W>(log(val.getStandard()), val.getDual()/val.getStandard());
}
template< typename W>
auto sqrt(const AutoDiff<W> &val){
    auto getsqrt=sqrt(val.getStandard());
    return AutoDiff<W>(getsqrt, val.getDual()/(2*getsqrt));
}
template<typename W>
auto erf(const AutoDiff<W> &val){
    return AutoDiff<W>(erf(val.getStandard()), (2.0/sqrt(M_PI))*exp(-val.getStandard()*val.getStandard())*val.getDual());
}
template<typename W>
auto abs(const AutoDiff<W> &val){
    return abs(val.getStandard());
}
template<typename W>
auto fabs(const AutoDiff<W> &val){
    return fabs(val.getStandard());
}
template<typename W,typename Z>
auto operator+(const AutoDiff<W> &left, const AutoDiff<Z> &right){
    return left.add(right);
}



template<typename W,typename Z>
auto operator+(const AutoDiff<W> &left, const Z& right){
    return left.add(right);
}
template< typename W,typename Z>
auto operator+(const W& left, const AutoDiff<Z> &right){
    return right.add(left);
}
template<typename W,typename Z>
auto operator-(const AutoDiff<W> &left, const AutoDiff<Z> &right){
    return left.subtract(right);
}
template<typename W>
auto operator-(const AutoDiff<W> &right){
    return AutoDiff<W>(-right.getStandard(), -right.getDual());
}
template<typename W,typename Z>
auto operator-(const AutoDiff<W> &left, const Z& right){
    return left.subtract(right);
}
template< typename W,typename Z>
auto operator-(const W& left, const AutoDiff<Z> &right){
    return right.multiply(-1).add(left);
}
template<typename W,typename Z>
auto operator*(const AutoDiff<W> &left, const AutoDiff<Z> &right){
    return left.multiply(right);
}
template<typename W,typename Z>
auto operator*(const AutoDiff<W> &left, const Z& right){
    return left.multiply(right);
}
template< typename W,typename Z>
auto operator*(const W& left, const AutoDiff<Z> &right){
    return right.multiply(left);
}
template<typename W,typename Z>
auto operator/(const AutoDiff<W> &left, const AutoDiff<Z> &right){
    return left.multiply(right.recipricol());
}
template<typename W,typename Z>
auto operator/(const AutoDiff<W> &left, const Z& right){
    return left.multiply(1.0/right);
}
template<typename W,typename Z>
auto operator/(const W& left, const AutoDiff<Z> &right){
    return right.recipricol().multiply(left);
}


template<typename U, typename W>
bool operator==(const AutoDiff<U> &left, const AutoDiff<W> &right){
    return left.getStandard()==right.getStandard() && left.getDual()==right.getDual();
}
template<typename U, typename W>
bool operator==(const U& left, const AutoDiff<W> &right){
    return left==right.getStandard();
}
template<typename U, typename W>
bool operator==(const AutoDiff<U> &left, const W& right){
    return left.getStandard()==right;
}
template<typename U, typename W>
bool operator!=(const AutoDiff<U> &left, const AutoDiff<W> &right){
    return left.getStandard()!=right.getStandard() && left.getDual()!=right.getDual();
}
template<typename U, typename W>
bool operator!=(const U& left, const AutoDiff<W> &right){
    return left!=right.getStandard();
}
template<typename U, typename W>
bool operator!=(const AutoDiff<U> &left, const W& right){
    return left.getStandard()!=right;
}
template<typename U, typename W>
bool operator>(const AutoDiff<U> &left, const AutoDiff<W> &right){
    return left.getStandard()>right.getStandard();
}
template<typename U, typename W>
bool operator>(const U& left, const AutoDiff<W> &right){
    return left>right.getStandard();
}
template<typename U, typename W>
bool operator>(const AutoDiff<U> &left, const W& right){
    return left.getStandard()>right;
}
template<typename U, typename W>
bool operator<(const AutoDiff<U> &left, const AutoDiff<W> &right){
    return left.getStandard()<right.getStandard();
}
template<typename U, typename W>
bool operator<(const U& left, const AutoDiff<W> &right){
    return left<right.getStandard();
}
template<typename U, typename W>
bool operator<(const AutoDiff<U> &left, const W& right){
    return left.getStandard()<right;
}

#endif