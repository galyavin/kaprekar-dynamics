// Independent whole finite graph and prefix checks by digit-column subtraction.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>
using Q=std::array<int,3>;
long long calls=0;
Q raw(int B,Q q){std::array<int,7>d={q[0],q[1],q[2],0,0,0,0},out{};int borrow=0;for(int j=6;j>=0;j--){int z=d[j]-d[6-j]-borrow;borrow=z<0;out[j]=z+borrow*B;}std::sort(out.begin(),out.end());calls++;return {out[6]-out[0],out[5]-out[1],out[4]-out[2]};}
long long index(Q q){long long u=q[0],v=q[1];return u*(u+1)*(u+2)/6+v*(v+1)/2+q[2]-1;}
int V(int B,Q q){return std::abs(3*B-4*q[0])+std::abs(4*(q[0]-q[1])-B)+std::abs(4*(q[1]-q[2])-B)+std::abs(4*q[2]-B);}
void require(bool b,const char*s){if(!b){std::cerr<<s<<'\n';std::exit(1);}}
int main(int argc,char**argv){int B=argc>1?std::atoi(argv[1]):400;bool prefixes=argc>2;long long size=(long long)B*(B+1)*(B+2)/6-1;std::vector<int> times(size,-1);int maxT=0,maxCycleV=0;long long cycleStates=0,positiveSources=0;Q witness{};
for(int u=1;u<B;u++)for(int v=0;v<=u;v++)for(int w=0;w<=v;w++){
 Q start={u,v,w};if(times[index(start)]>=0)continue;
 std::vector<Q> path;Q q=start;
 while(times[index(q)]==-1){times[index(q)]=-2;path.push_back(q);q=raw(B,q);require(q[0]>0,"Constant image");}
 if(times[index(q)]==-2){auto it=std::find(path.begin(),path.end(),q);require(it!=path.end(),"Cycle missing in own path");for(auto jt=it;jt!=path.end();jt++){times[index(*jt)]=0;cycleStates++;maxCycleV=std::max(maxCycleV,V(B,*jt));require((*jt)[2]>0,"Zero cycle");}path.erase(it,path.end());}
 int t=times[index(q)];for(auto it=path.rbegin();it!=path.rend();it++){times[index(*it)]=++t;if(t>maxT){maxT=t;witness=*it;}}
}
require(maxCycleV<=272,"Cycle potential bound");require(maxT<=54LL*B+614131,"Global pair clock bound");
if(prefixes)for(int u=1;u<B;u++)for(int v=0;v<=u;v++)for(int w=0;w<=v;w++){
 Q start={u,v,w},q=start;int p0=V(B,start);bool positive=w>0;
 for(int j=1;j<=9;j++){Q next=raw(B,q);if(q[2]>0)require(V(B,next)<=V(B,q)+8,"One-step increase");q=next;positive=positive&&q[2]>0;if(j>=6)require(q[2]>0,"Six-step permanent positivity prefix");}
 if(positive){positiveSources++;if(p0>200)require(V(B,q)<p0,"Nine-step strict decrease");}
}
std::cout<<"{\"base\":"<<B<<",\"states\":"<<size<<",\"raw_column_updates\":"<<calls<<",\"cycle_states\":"<<cycleStates<<",\"maximum_pair_time\":"<<maxT<<",\"maximum_cycle_potential\":"<<maxCycleV<<",\"maximum_time_source\":["<<witness[0]<<','<<witness[1]<<','<<witness[2]<<"],\"positive_nine_prefixes\":"<<positiveSources<<",\"all_prefixes_checked\":"<<(prefixes?"true":"false")<<"}\n";
}
