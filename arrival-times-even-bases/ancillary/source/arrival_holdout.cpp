// Independent full pair-graph check using sorted subtraction output digits.
// Usage: arrival_holdout core exponent fixed
// No inverse branches, symbolic boundary rules, or limiting formulas are used.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using U = std::uint64_t;
using Pair = std::pair<int,int>;
U index_of(int u, int v) { return U(u)*(u+1)/2+v-1; }

Pair step(int b, Pair p) {
    int u=p.first, v=p.second;
    std::array<int,4> ds;
    if (v) ds={u,v-1,b-1-v,b-u};
    else ds={u-1,b-1,b-1,b-u};
    std::sort(ds.begin(),ds.end());
    return {ds[3]-ds[0],ds[2]-ds[1]};
}

int main(int argc, char** argv) {
    if (argc!=4) throw std::runtime_error("Expected core, exponent, fixed");
    const int m=std::stoi(argv[1]), n=std::stoi(argv[2]);
    if (n<0 || n>15 || m<3 || m%2==0) throw std::runtime_error("Invalid parameters");
    const int scale=1<<n, b=m*scale;
    if(b>6000) throw std::runtime_error("Holdout size limited to base 6000");
    const std::string kind=argv[3];
    std::set<Pair> targets;
    if(kind!="fixed" || m!=5) throw std::runtime_error("Expected the base-five fixed point");
    targets.insert({3*scale,scale});
    for(auto p:targets) if(!targets.count(step(b,p))) throw std::runtime_error("Invalid target");
    const U count=U(b)*(b+1)/2-1;
    // -2: unknown; -3: current path; -1: misses the target.
    std::vector<int> distance(count,-2);
    for(auto p:targets) distance[index_of(p.first,p.second)]=0;
    std::vector<U> path;
    std::map<int,U> histogram;
    U basin_pairs=0, all_weight=0;
    for(int u=1;u<b;++u) for(int v=0;v<=u;++v) {
        Pair q={u,v};
        path.clear();
        while(distance[index_of(q.first,q.second)]==-2) {
            U idx=index_of(q.first,q.second);
            distance[idx]=-3; path.push_back(idx);
            q=step(b,q);
        }
        int value=distance[index_of(q.first,q.second)];
        if(value<0) value=-1;
        for(auto it=path.rbegin();it!=path.rend();++it) {
            if(value>=0) ++value;
            distance[*it]=value;
        }
        U weight=!v ? U(b-u)*(12*u-4) :
            (u==v ? U(6)*(b-u) : U(24)*(b-u)*(u-v));
        all_weight+=weight;
        int t=distance[index_of(u,v)];
        if(t>=0) { histogram[t+1]+=weight; ++basin_pairs; }
    }
    if(!targets.empty()) { histogram[1]-=targets.size(); histogram[0]+=targets.size(); }
    const U eligible=U(b)*b*b*b-b;
    if(all_weight!=eligible) throw std::runtime_error("Input weights do not partition inputs");
    U total=0, first=0, second=0;
    for(auto [t,w]:histogram) { total+=w; first+=t*w; second+=U(t)*t*w; }
    std::cout<<"{\"core\":"<<m<<",\"exponent\":"<<n<<",\"base\":"<<b
             <<",\"checked_pairs\":"<<count<<",\"basin_pairs\":"<<basin_pairs
             <<",\"basin_input_count\":"<<total<<",\"complement_input_count\":"<<eligible-total
             <<",\"first_moment_sum\":"<<first<<",\"second_moment_sum\":"<<second
             <<",\"input_arrival_counts\":{";
    bool comma=false;
    for(auto [t,w]:histogram) if(w) {
        if(comma) std::cout<<",";
        comma=true; std::cout<<"\""<<t<<"\":"<<w;
    }
    std::cout<<"}}\n";
}
