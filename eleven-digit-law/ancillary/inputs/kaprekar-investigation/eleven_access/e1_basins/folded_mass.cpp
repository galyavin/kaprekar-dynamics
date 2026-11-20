// Rigorous lower weights on the original leading population, using boxes of
// folded lows after the initial sort. Arithmetic is entirely integral.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using Wide = __int128_t;
using Vec = std::array<int,5>;
using Matrix = std::array<Vec,5>;

struct Form {
    int center;
    Vec coefficients;
    int source;
};

Vec add(const Vec& a, const Vec& b, int sign = 1) {
    Vec result{};
    for (int i=0;i<5;++i) result[i] = a[i] + sign*b[i];
    return result;
}

int radius(const Vec& a) {
    int total=0;
    for (int c : a) total += std::abs(c);
    return total;
}

int first_entry(const Vec& bins, int m, int through) {
    const int unit=4*m;
    Vec centers{};
    Matrix basis{}, gradients{};
    for (int j=0;j<5;++j) basis[j][j]=1;
    for (int j=0;j<5;++j) {
        centers[j]=unit-(2*bins[j]+1)-(j ? 2*bins[j-1]+1 : 0);
        gradients[j]=add(Vec{},basis[j],-1);
        if (j) gradients[j]=add(gradients[j],basis[j-1],-1);
    }
    const Vec positive{1,-5,2,-4,3};
    int consecutive=0;
    // The initial fold produces the time-one difference vector.
    for (int t=1;t<=through+4;++t) {
        std::array<Form,5> lows{};
        for (int j=0;j<5;++j) {
            const int sign = 2*centers[j]>unit ? 1 : -1;
            if (std::abs(2*centers[j]-unit) < 2*radius(gradients[j])) return -1;
            lows[j]={sign>0 ? unit-centers[j] : centers[j],
                     add(Vec{},gradients[j],-sign),(j+1)*sign};
        }
        std::sort(lows.begin(),lows.end(),[](const Form& x,const Form& y) {
            return x.center < y.center;
        });
        bool positive_order=true;
        for(int j=0;j<5;++j) {
            if (j && lows[j].center-lows[j-1].center
                       < radius(add(lows[j].coefficients,lows[j-1].coefficients,-1)))
                return -1;
            positive_order = positive_order && lows[j].source == positive[j];
        }
        consecutive = positive_order ? consecutive+1 : 0;
        if (consecutive == 5) return t-4;
        for (int j=0;j<5;++j) {
            centers[j]=unit-lows[j].center-(j ? lows[j-1].center : 0);
            gradients[j]=add(Vec{},lows[j].coefficients,-1);
            if(j) gradients[j]=add(gradients[j],lows[j-1].coefficients,-1);
        }
    }
    return 0;
}

Wide integrated_numerator(const Vec& bins, int m) {
    // On this box X_i=24*m*low_i has integer mean and variance 12.
    // Along each side, d=E[gap product], e=E[last coordinate*gap product].
    // The two sides use disjoint variables, so the middle gap is exact.
    Wide result=0;
    for(int mask=0;mask<16;++mask) {
        std::array<Wide,2> d{1,1},e{0,0};
        for(int j=0;j<5;++j) {
            const int side=(mask>>j)&1;
            const Wide center=6*(2*bins[j]+1);
            const Wide next_e=(center*center+12)*d[side]-center*e[side];
            d[side]=center*d[side]-e[side];
            e[side]=next_e;
        }
        const Wide mass=24*m*d[0]*d[1]-e[0]*d[1]-d[0]*e[1];
        if (mass<=0) std::abort();
        result += 2*mass;  // Reflection swaps the two sides.
    }
    return result;
}

void check_integration(const Vec& bins,int m) {
    // Tensor Simpson is exact independently, since each coordinate occurs
    // in at most two gap factors. Here all original masks are evaluated.
    const std::array<int,3> weights{1,4,1};
    Wide weighted=0;
    for(int node=0;node<243;++node) {
        Vec lows{};
        int digits=node;
        int weight=1;
        for(int j=0;j<5;++j) {
            const int k=digits%3;
            digits/=3;
            lows[j]=6*(2*bins[j]+k);
            weight*=weights[k];
        }
        for(int mask=0;mask<32;++mask) {
            Vec q{};
            for(int j=0;j<5;++j)
                q[j] = (mask>>j)&1 ? 24*m-lows[j] : lows[j];
            std::sort(q.begin(),q.end(),std::greater<int>());
            Wide product=q.back();
            int previous=24*m;
            for(int x:q) {
                product*=previous-x;
                previous=x;
            }
            weighted+=weight*product;
        }
    }
    if(weighted != 7776*integrated_numerator(bins,m)) std::abort();
}

std::string decimal(Wide n) {
    if (!n) return "0";
    std::string digits;
    while(n) {
        digits.push_back(char('0'+n%10));
        n/=10;
    }
    std::reverse(digits.begin(),digits.end());
    return digits;
}

int main(int argc, char** argv) {
    if (argc!=3) return 1;
    const int m=std::atoi(argv[1]), through=std::atoi(argv[2]);
    // With through<=19, worst-case affine coefficients and integer centres fit
    // int. With m<=256 the six-factor integrals and their sums fit Wide.
    if(m<5 || m>256 || through<5 || through>19) return 1;
    check_integration({0,1,2,3,4},m);
    check_integration({m-5,m-4,m-3,m-2,m-1},m);
    if(m>=16) check_integration({0,m/4,m/2,3*m/4,m-1},m);
    std::vector<Wide> totals(through+1);
    std::vector<std::uint64_t> counts(through+1);
    std::uint64_t ambiguous=0, outside=0;
    for(int a=0;a<m;++a)
        for(int b=a+1;b<m;++b)
            for(int c=b+1;c<m;++c)
                for(int d=c+1;d<m;++d)
                    for(int e=d+1;e<m;++e) {
                        const Vec bins{a,b,c,d,e};
                        const int t=first_entry(bins,m,through);
                        if(t<0) ++ambiguous;
                        else if(!t) ++outside;
                        else {
                            ++counts[t];
                            totals[t]+=integrated_numerator(bins,m);
                        }
                    }
    std::cout << "grid " << m << " through " << through << " ambiguous " << ambiguous
              << " outside " << outside << '\n';
    for(int t=1;t<=through;++t)
        std::cout << t << ' ' << counts[t] << ' ' << decimal(totals[t]) << '\n';
}
