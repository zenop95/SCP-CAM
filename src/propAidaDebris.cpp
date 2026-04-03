#include <nlohmann/json.hpp>
#include <vector>
#include <array>
#include <fstream>
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"
#include "dynorb/AIDA.h"
#include "dynorb/AIDAwrappers.h"
#include "prop_utils.h"
using json = nlohmann::json;
using namespace DACE;
using namespace props::io;

std::vector<std::vector<double>>
to_vec2m(const AlgebraicMatrix<double>& M, int R, int C) {
    std::vector<std::vector<double>> out(R, std::vector<double>(C));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
            out[r][c] = M.at(r, c);
    return out;
}

int main(){
 try{
  long t0=now_ms(); long time_acc_ms=0;
  auto jin=read_json("./data_sharing/initial_state.json");
  int N=jin["N"].get<int>();
  double Lsc=jin["scaling"]; double et=jin["et"];
  std::vector<double> t=jin["t"].get<std::vector<double>>();
  auto x0v=jin["x0"].get<std::vector<double>>();

  DA::init(1,6); 
  DA::setEps(1e-15);

  const auto aida = read_aida_params_json("./data_sharing/AIDA_init.json");
  load_spice_kernel("./data/kernel.txt");
  auto coeff=compute_aida_coeff(aida);
  int flags[3]={aida.flag1,aida.flag2,aida.flag3};
  AIDAScaledDynamics<DA> dyn(
        "./data/gravmodels/egm2008",
        aida.gravOrd,
        flags,
        coeff.Bfactor,
        coeff.SRPC
    );
  bool isNotKep = (aida.gravOrd > 1) || (aida.flag1 > 0) || (aida.flag2 > 0) || (aida.flag3 > 0);

  AlgebraicVector<DA> x0(6),x(6),u(3);
  u={0.0,0.0,0.0}; for(int j=0;j<6;j++)x0[j]=x0v[j]+DA(j+1);

  // first node
  x0=RK78Sc(6,x0,u,et,et+t[0],1.0,Lsc,0,isNotKep,dyn);
  for(int j=0;j<6;j++)x0[j]=cons(x0[j])+DA(j+1);
  AlgebraicMatrix<double> STM=astro::stmDace(x0,6,6);

  std::vector<std::array<double,6>> constPart; std::vector<std::vector<std::vector<double>>> maps;
  {
    std::array<double,6> row{}; for(int j=0;j<6;j++)row[j]=cons(x0[j]);
    constPart.push_back(row); maps.push_back(to_vec2m(STM,6,6));
  }

  long t1=now_ms(); time_acc_ms+=now_ms()-t1;

  // remaining nodes
  for(int i=1;i<N;i++){
    double t0n=et+t[i-1], tfn=et+t[i];
    x=RK78Sc(6,x0,u,t0n,tfn,1.0,Lsc,0,isNotKep,dyn);
    for(int j=0;j<6;j++){
        x0[j] = cons(x[j]) + DA(j+1);
        STM   = astro::stmDace(x,6,6);
    }

    std::array<double,6> row{}; for(int j=0;j<6;j++)row[j]=cons(x[j]);
    constPart.push_back(row); maps.push_back(to_vec2m(STM,6,6));
  }

  time_acc_ms+=now_ms()-t0;
  json jout; jout["constPart"]=constPart; jout["stm"]=maps; jout["timeMs"]=time_acc_ms;
  write_json("./data_sharing/out_prop_sec.json",jout);
  clear_spice_kernels(); return 0;
 } catch(const std::exception& ex){fprintf(stderr,"Error: %s\n",ex.what()); return 1;}
}