// Direct column subtraction from one original inverse branch of the separating route.
// No core, affine-return, or certificate code is used here.
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using State = std::array<int,4>;
// Rejection on nine uniform bits makes the seeded input sequence portable.
int draw_offset(std::mt19937& rng) {
    unsigned int bits;
    do {
        bits = rng() & 511u;
    } while (bits >= 421);
    return static_cast<int>(bits) - 210;
}

State step(int base, const State& q) {
    std::array<int,9> digits{};
    int borrow=0;
    for(int i=0;i<9;++i) {
        int column=(i<4 ? -q[i] : i==4 ? 0 : q[8-i])-borrow;
        borrow=column<0;
        digits[i]=column+base*borrow;
    }
    std::sort(digits.begin(), digits.end());
    State out;
    for(int i=0;i<4;++i) out[i]=digits[8-i]-digits[i];
    return out;
}
int main() {
    std::mt19937 rng(1011711);
    const double weights[5][4]={
        {11./48,5./48,137./288,55./288},{29./128,5./48,47./96,23./128},
        {29./144,5./48,607./1152,193./1152},{107./576,5./48,69./128,197./1152},
        {233./1152,5./48,73./144,215./1152}};
    const double rates[]={3.5,4.5,50./9,37./6};
    const int samples=4000;
    for(int residue=0;residue<5;++residue) {
        int base=300035+residue;
        std::array<int,4> counts{};
        double largest_clock_error=0, largest_terminal_error=0;
        for(int i=0;i<samples;++i) {
            State q;
            const int numer[]={43,39,31,28};
            for(int j=0;j<4;++j) q[j]=(base*numer[j])/50+draw_offset(rng);
            // All four x_i exceed 1/2 and decrease strictly here.
            double u=double(q[0])/base, v=double(q[0]+q[1])/base-1,
                   z=double(q[2]+q[3])/base-1;
            double L=1.75*u-v-.25*z-.75, s=-.25*u+.5*v-.25*z+.75;
            // Floyd detection plus a second pass recovers the first cyclic pair.
            State slow=step(base,q),fast=step(base,step(base,q));
            while(slow!=fast){slow=step(base,slow);fast=step(base,step(base,fast));}
            int h=0;
            for(slow=q;slow!=fast;++h){slow=step(base,slow);fast=step(base,fast);}
            if(h==0) throw std::runtime_error("unexpected cyclic source");
            State terminal=step(base,slow);  // differences of the first cyclic number
            double time=double(h+1)/base, best=1e30, second=1e30, terr=0, zerr=0;
            int selected=-1;
            for(int j=0;j<4;++j) {
                double a=j==0?(s+4)/6:.8;
                std::array<double,4> endpoint={a,2*a-1,3*a-2,1-a};
                double dt=std::abs(time-(L+rates[j]*(s-.8))),dz=0;
                for(int k=0;k<4;++k) dz=std::max(dz,std::abs(double(terminal[k])/base-endpoint[k]));
                double error=std::max(dt,dz);
                if(error<best){second=best;best=error;selected=j;terr=dt;zerr=dz;}
                else second=std::min(second,error);
            }
            if(best>.0015 || second<2*best) throw std::runtime_error("ambiguous or inaccurate atom");
            ++counts[selected];
            largest_clock_error=std::max(largest_clock_error,terr);
            largest_terminal_error=std::max(largest_terminal_error,zerr);
        }
        std::cout<<base;
        for(int j=0;j<4;++j) {
            std::cout<<" "<<counts[j];
            if(std::abs(double(counts[j])/samples-weights[residue][j])>.035)
                throw std::runtime_error("large frequency discrepancy");
        }
        std::cout<<" max clock and terminal errors "<<largest_clock_error<<" "<<largest_terminal_error<<std::endl;
    }
}
