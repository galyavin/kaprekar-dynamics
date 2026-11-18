// Independent finite checks using nine digit-column subtraction.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
using Q=std::array<int,4>;
long long calls=0;
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
Q raw(int B,Q q){
 std::array<int,9>d={q[0],q[1],q[2],q[3],0,0,0,0,0},out{};
 int borrow=0;
 for(int j=8;j>=0;j--){int z=d[j]-d[8-j]-borrow;borrow=z<0;out[j]=z+borrow*B;}
 require(borrow==0,"Final borrow");std::sort(out.begin(),out.end());calls++;
 return {out[8]-out[0],out[7]-out[1],out[6]-out[2],out[5]-out[3]};
}
std::array<int,5> gaps(int B,Q q){return {B-q[0],q[0]-q[1],q[1]-q[2],q[2]-q[3],q[3]};}
int G(int B,Q q){auto g=gaps(B,q);int s=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)s+=std::abs(g[i]-g[j]);return s;}
int P(int B,Q q){auto g=gaps(B,q);int m=*std::min_element(g.begin(),g.end()),M=*std::max_element(g.begin(),g.end());return B-4*m-std::min(M,2*m);}
bool cover(int B,Q q){
 int u=q[0],v=q[1],w=q[2];
 if(!v)return true;
 if(!w)return std::abs(u-v)<=3||std::abs(u+v-B)<=3||std::abs(2*u-B)<=4||std::abs(2*v-B)<=4||B-u<=3;
 if(v==w&&(std::abs(u-v)<=2||std::abs(u+v-B)<=2||std::abs(2*v-B)<=3))return true;
 if(u==v&&(std::abs(u+w-B)<=2||std::abs(2*u-B)<=3))return true;
 if(std::abs(2*u-B)<=3&&std::abs(2*v-B)<=3)return true;
 if(B-u<=3&&(std::abs(2*v-B)<=3||v==w||std::abs(v+w-B)<=3))return true;
 if(std::abs(4*u-3*B)<=8)for(int x:{v,w})for(int j:{1,2,3})if(std::abs(4*x-j*B)<=12)return true;
 return false;
}
long long index(Q q){long long u=q[0],v=q[1],w=q[2];return u*(u+1)*(u+2)*(u+3)/24+v*(v+1)*(v+2)/6+w*(w+1)/2+q[3]-1;}
void zero_checks(int B){
 long long sources=0,covered=0,returning=0,early=0;int longest=0,maxreturn=0;
 for(int u=1;u<B;u++)for(int v=0;v<=u;v++)for(int w=0;w<=v;w++){
  Q start={u,v,w,0},q=start;bool membership=cover(B,start);covered+=membership;sources++;
  for(int j=1;j<=4;j++){q=raw(B,q);if(!q[3]){returning++;maxreturn=std::max(maxreturn,j);require(membership,"Uncovered return source");break;}}
  q=start;std::vector<Q> path;int positive=0,clock=0;
  while(true){
   auto it=std::find(path.begin(),path.end(),q);
   if(it!=path.end()){early++;clock=it-path.begin();break;}
   path.push_back(q);positive=q[3]?positive+1:0;
   if(positive==4){clock=path.size()-4;break;}
   q=raw(B,q);require(path.size()<10000,"Finite-check safety cap");
  }
  longest=std::max(longest,clock);
 }
 require(covered<=308LL*B,"Cover count bound");
 require(longest<=1232LL*B+7,"Zero preparation bound");
 std::cout<<"{\"mode\":\"zero\",\"base\":"<<B<<",\"sources\":"<<sources<<",\"covered\":"<<covered<<",\"returning\":"<<returning<<",\"largest_first_return\":"<<maxreturn<<",\"early_cycles\":"<<early<<",\"maximum_clearance_or_capture\":"<<longest<<",\"updates\":"<<calls<<"}\n";
}
void graph_checks(int B){
 long long size=(long long)B*(B+1)*(B+2)*(B+3)/24-1,positive=0,tested=0;
 std::vector<int> times(size,-1);int maxT=0,maxcycleP=0,maxincrease=-100,latestdescent=0;long long cycles=0;
 for(int u=1;u<B;u++)for(int v=0;v<=u;v++)for(int w=0;w<=v;w++)for(int z=0;z<=w;z++){
  Q q={u,v,w,z};if(times[index(q)]>=0)continue;std::vector<Q> path;
  while(times[index(q)]==-1){times[index(q)]=-2;path.push_back(q);q=raw(B,q);}
  if(times[index(q)]==-2){auto it=std::find(path.begin(),path.end(),q);for(auto jt=it;jt!=path.end();jt++){times[index(*jt)]=0;cycles++;if((*jt)[3])maxcycleP=std::max(maxcycleP,P(B,*jt));}path.erase(it,path.end());}
  int t=times[index(q)];for(auto it=path.rbegin();it!=path.rend();it++){times[index(*it)]=++t;maxT=std::max(maxT,t);}
 }
 for(int u=1;u<B;u++)for(int v=1;v<=u;v++)for(int w=1;w<=v;w++)for(int z=1;z<=w;z++){
  Q start={u,v,w,z};int initial=G(B,start);Q q=raw(B,start);positive++;
  int increase=G(B,q)-initial;maxincrease=std::max(maxincrease,increase);require(increase<=8,"One-step dispersion increase");
  if(P(B,start)<=100)continue;
  bool success=false,zero=false;
  for(int j=1;j<=9;j++){
   if(!q[3]){zero=true;break;}
   if(G(B,q)<initial){latestdescent=std::max(latestdescent,j);success=true;break;}
   if(j<9)q=raw(B,q);
  }
  if(!zero){tested++;require(success,"Nine-step lower-hit candidate fails");}
 }
 std::cout<<"{\"mode\":\"graph\",\"base\":"<<B<<",\"states\":"<<size<<",\"positive_sources\":"<<positive<<",\"candidate_descent_sources\":"<<tested<<",\"latest_descent\":"<<latestdescent<<",\"maximum_one_step_increase\":"<<maxincrease<<",\"maximum_pair_clock\":"<<maxT<<",\"cycle_states\":"<<cycles<<",\"maximum_positive_cycle_residual\":"<<maxcycleP<<",\"updates\":"<<calls<<"}\n";
}
int main(int argc,char**argv){try{require(argc==3,"Usage: binary base zero|graph");int B=std::atoi(argv[1]);require(B>=40,"B>=40 required");if(std::string(argv[2])=="zero")zero_checks(B);else{require(std::string(argv[2])=="graph","Unknown mode");graph_checks(B);}}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
