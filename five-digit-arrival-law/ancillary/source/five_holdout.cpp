// Independent five-slot reference. No borrowed-digit branches or excursions.
// Build: c++ -O2 -std=c++17 five_holdout.cpp -o /tmp/five-holdout
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using Pair = std::pair<int,int>;
using Wide = unsigned __int128;

std::string decimal(Wide n) {
    if (!n) return "0";
    std::string s;
    while (n) { s.push_back(char('0' + n % 10)); n /= 10; }
    std::reverse(s.begin(), s.end()); return s;
}

Pair step(int base, Pair p) {
    // (u,v,0,0,0) is a sorted representative with differences (u,v).
    std::array<int,5> d{p.first,p.second,0,0,0};
    std::uint64_t down=0,up=0;
    for (int i=0;i<5;i++) { down=down*base+d[i]; up=up*base+d[4-i]; }
    auto n=down-up;
    for (int i=0;i<5;i++) { d[i]=n%base; n/=base; }
    std::sort(d.begin(),d.end());
    return {d[4]-d[0],d[3]-d[1]};
}

Wide weight(int base, Pair p) {
    auto [u,v]=p;
    // Count c at an endpoint or inside [d,b], and d at an endpoint or
    // inside [e,a-v]. These multinomial categories are independently
    // checked by actual ordered-input enumeration in five_checks.py.
    Wide per_shift;
    if (!v) per_shift=20*(u-1)+5+5;
    else if (u==v) per_shift=30*(v-1)+10+10;
    else per_shift=Wide(u-v-1)*(120*(v-1)+60+60)
                       +2*Wide(60*(v-1)+20+30);
    return (base-u)*per_shift;
}

int main(int argc,char** argv) {
    if (argc!=2) { std::cerr << "usage: five-holdout BASE\n"; return 2; }
    int b=std::stoi(argv[1]);
    if (b<2 || b>2000) throw std::invalid_argument("2 <= base <= 2000");
    auto id=[](Pair p){return p.first*(p.first+1)/2+p.second-1;};
    int n=(b-1)*(b+2)/2;
    std::vector<Pair> pairs(n);
    for(int u=1;u<b;u++) for(int v=0;v<=u;v++) pairs[id({u,v})]={u,v};
    std::vector<int> next(n),dest(n,-1),distance(n),mark(n,-1),position(n);
    for(int i=0;i<n;i++) next[i]=id(step(b,pairs[i]));
    std::vector<std::vector<int>> cycles;
    for(int start=0;start<n;start++) {
        if(dest[start]>=0)continue;
        std::vector<int> path;int q=start;
        while(dest[q]<0 && mark[q]!=start) {
            mark[q]=start;position[q]=path.size();path.push_back(q);q=next[q];
        }
        if(dest[q]<0) {
            int cut=position[q], c=cycles.size();
            std::vector<int> cycle(path.begin()+cut,path.end());
            auto smallest=std::min_element(cycle.begin(),cycle.end());
            std::rotate(cycle.begin(),smallest,cycle.end());
            cycles.push_back(cycle);
            for(int p:cycle){dest[p]=c;distance[p]=0;}
            path.resize(cut);
        }
        for(auto it=path.rbegin();it!=path.rend();it++) {
            dest[*it]=dest[next[*it]];distance[*it]=distance[next[*it]]+1;
        }
    }
    std::vector<std::map<int,Wide>> hist(cycles.size());
    Wide total=0,rectangle=0,tail=0,first=0,second=0;
    for(int i=0;i<n;i++) {
        auto w=weight(b,pairs[i]);hist[dest[i]][distance[i]+1]+=w;total+=w;
        auto [u,v]=pairs[i];
        if(5*u>=4*b && 20*u<=17*b && 5*v>=2*b && 20*v<=9*b)rectangle+=w;
    }
    for(unsigned c=0;c<cycles.size();c++) {
        hist[c][1]-=cycles[c].size();hist[c][0]+=cycles[c].size();
        for(auto [t,w]:hist[c]) {
            first+=t*w;second+=Wide(t)*t*w;if(25*t>=b)tail+=w;
        }
    }
    std::cout << "{\"base\":"<<b<<",\"pair_count\":"<<n
      <<",\"total_input_weight\":"<<decimal(total)
      <<",\"rectangle_weight\":"<<decimal(rectangle)
      <<",\"tail_at_B_over_25\":"<<decimal(tail)
      <<",\"first_moment_sum\":"<<decimal(first)
      <<",\"second_moment_sum\":"<<decimal(second)<<",\"cycles\":[";
    for(unsigned c=0;c<cycles.size();c++) {
        if(c)std::cout<<',';std::cout<<"{\"pairs\":[";
        for(unsigned j=0;j<cycles[c].size();j++){
            if(j)std::cout<<',';auto p=pairs[cycles[c][j]];
            std::cout<<'['<<p.first<<','<<p.second<<']';
        }
        std::cout<<"],\"input_arrival_counts\":{";bool comma=false;
        for(auto [t,w]:hist[c])if(w){
            if(comma)std::cout<<',';comma=true;
            std::cout<<'"'<<t<<"\":"<<decimal(w);
        }
        std::cout<<"}}";
    }
    std::cout<<"]}\n";
}
