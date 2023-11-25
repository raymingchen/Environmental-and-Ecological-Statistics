from numpy import *
from math import *  ##for factorial
from tkinter import filedialog as fd
import pandas as pd
from itertools import *  ## for permutations
import array_to_latex as atl ##for output presentations
from matplotlib import pyplot
from tabulate import tabulate
## processing preliminary input data
csv_data = pd.read_csv(fd.askopenfilename())
InputData_opp=array(csv_data)
csv_data = pd.read_csv(fd.askopenfilename())
InputData_risk=array(csv_data)
##Mm2 is a nature-numbered m-by-2 matrix, in which m is the number of
##  inputs. Two columns correspond to two variables
## p is a natural number (including 0)
def PtoPposCauslity(Mm22,p): ## point to point positive causality
  indepenVar=Mm22[:,0];
  depenVar=Mm22[:,1];
  numOfCorrEles=sum(indepenVar==depenVar)
  indepenVar=indepenVar.tolist()
  depenVar = depenVar.tolist()
  p1 = indepenVar.index(p)
  p2 = depenVar.index(p)
  if p1!=p2:
    ptoPposCauslity=sign(p2-p1)/(2**abs(p2-p1))
  else:
      if numOfCorrEles<=1:ptoPposCauslity=0
      else: ptoPposCauslity=1
  return ptoPposCauslity

##
def PtoPnegCauslity(Mm22,p): ## point to point negative causality
  indepenVar=Mm22[:,0];
  depenVar=Mm22[:,1];
  numOfCorrEles=sum(indepenVar==(shape(Mm22)[0]-1-depenVar))
  indepenVar=indepenVar.tolist()
  depenVar = depenVar.tolist()
  p1 = indepenVar.index(p)
  p2 = depenVar.index(shape(Mm22)[0]-1-p)
  if p1!=p2:
    ptoPnegCauslity=-sign(p2-p1)/(2**abs(p2-p1))
  else:
      if numOfCorrEles<=1:ptoPnegCauslity=0
      else: ptoPnegCauslity=-1
  return ptoPnegCauslity
##Mm2::any real-numbered m-by-2 matrix
def TwoColposCausality(Mm2):
  Mm21=Mm2.T.argsort().argsort(); Mm21=Mm21.T
  ind=Mm21[:,0].argsort()
  Mm22=Mm21[ind,:];Mm22=array(Mm22,dtype=int)
  m=shape(Mm22)[0]
  vToVposCausality=ones(m)
  for i in range(m):
    vToVposCausality[i]=PtoPposCauslity(Mm22,i)
    twoColposCausality=sum(vToVposCausality)/m
  return twoColposCausality

def TwoColnegCausality(Mm2):
  Mm21=Mm2.T.argsort().argsort(); Mm21=Mm21.T
  ind=Mm21[:,0].argsort()
  Mm22=Mm21[ind,:];Mm22=array(Mm22,dtype=int)
  m=shape(Mm22)[0]
  vToVposCausality=ones(m)
  for i in range(m):
    vToVposCausality[i]=PtoPnegCauslity(Mm22,i)
    twoColnegCausality=sum(vToVposCausality)/m
  return twoColnegCausality
## Co-causal matrix
## Mat is a m-by-n matrix, where n indicates n variables(colulm)
def CoposCausMat(Mat):
  m=shape(Mat)[0];n=shape(Mat)[1]
  coposcausmat=ones((n,n))
  for i in range(n):
    for j in range(n):
      Mm2=Mat[:,(i,j)]
      coposcausmat[i,j]=TwoColposCausality(Mm2)
  return(coposcausmat)
def ConegCausMat(Mat):
  m=shape(Mat)[0];n=shape(Mat)[1]
  conegcausmat=ones((n,n))
  for i in range(n):
    for j in range(n):
      Mm2=Mat[:,(i,j)]
      conegcausmat[i,j]=TwoColnegCausality(Mm2)
  return(conegcausmat)

def PosPDF(m):
  perm=[ele for ele in permutations(range(m+1))]
  posPdf=ones(factorial(m+1))
  i=0
  for ele in perm:
    Mm2=column_stack((range(m+1),ele))
    posPdf[i]=TwoColposCausality(Mm2)
    i=i+1
  return posPdf
def NegPDF(m):
  perm=[ele for ele in permutations(range(m+1))]
  negPdf=ones(factorial(m+1))
  i=0
  for ele in perm:
    Mm2=column_stack((range(m+1),ele))
    negPdf[i]=TwoColnegCausality(Mm2)
    i=i+1
  return negPdf
## examples
Mm2=array([[2.3,5],[3.6,1.7],[1.3,4.5],[5.9,3],[4.4,2],[5,0.5]])
res2=TwoColposCausality(Mm2)
print(res2)
##
Mm22=array([[0,5],[1,1],[2,2],[3,4],[4,3],[5,0]])
resu22_0=PtoPposCauslity(Mm22,0)
resu22_1=PtoPposCauslity(Mm22,1)
resu22_2=PtoPposCauslity(Mm22,2)
resu22_3=PtoPposCauslity(Mm22,3)
resu22_4=PtoPposCauslity(Mm22,4)
resu22_5=PtoPposCauslity(Mm22,5)
print('resu22(0),resu22(1),resu22(2),resu22(3),resu22(4),\n resu22(5)=',
      resu22_0,resu22_1,resu22_2,resu22_3,resu22_4,resu22_5)
##
Mm22=array([[0,5],[1,1],[2,2],[3,4],[4,3],[5,0]])
resu22_0=PtoPnegCauslity(Mm22,0)
resu22_1=PtoPnegCauslity(Mm22,1)
resu22_2=PtoPnegCauslity(Mm22,2)
resu22_3=PtoPnegCauslity(Mm22,3)
resu22_4=PtoPnegCauslity(Mm22,4)
resu22_5=PtoPnegCauslity(Mm22,5)
print('resu22(0),resu22(1),resu22(2),resu22(3),resu22(4),\n resu22(5)=',
      resu22_0,resu22_1,resu22_2,resu22_3,resu22_4,resu22_5)
##
Mm2=array([[2.3,5],[3.6,1.7],[1.3,4.5],[5.9,3],[4.4,2],[5,0.5]])
resul2=TwoColnegCausality(Mm2)
print(resul2)
##
Mm2_opp=InputData_opp
Mm2_risk=InputData_risk
res_opp=CoposCausMat(Mm2_opp)
res_risk=CoposCausMat(Mm2_risk)
print('res_opp',res_opp)
print('res_risk',res_risk)
##
Mm2_opp=InputData_opp
Mm2_risk=InputData_risk
res_opp=ConegCausMat(Mm2_opp)
res_risk=ConegCausMat(Mm2_risk)
print('res_opp',res_opp)
print('res_risk',res_risk)

## to latex
## res_opp_latex=atl.to_ltx(res_opp)
## df = pd.DataFrame(res_opp)
##
posPDF_6=PosPDF(6)
pyplot.hist(posPDF_6)
q05_pos_6=quantile(posPDF_6,0.05)
q95_pos_6=quantile(posPDF_6,0.95)
print('5% quantile of posPDF_6=',q05_pos_6)
print('95% quantile of posPDF_6=',q95_pos_6)
##
negPDF_6=NegPDF(6)
pyplot.hist(negPDF_6)
q05_neg_6=quantile(negPDF_6,0.05)
q95_neg_6=quantile(negPDF_6,0.95)
print('5% quantile of negPDF_6=',q05_neg_6)
print('95% quantile of negPDF_6=',q95_neg_6)

fig, axes = pyplot.subplots(1, 2)
fig.suptitle("Positive (left) and negative (right) PDF(6)")
axes[0].hist(posPDF_6)
axes[1].hist(negPDF_6)
##
posPDF_7=PosPDF(7)
pyplot.hist(posPDF_7)
q05_pos_7=quantile(posPDF_7,0.05)
q95_pos_7=quantile(posPDF_7,0.95)
print('5% quantile of posPDF_7=',q05_pos_7)
print('95% quantile of posPDF_7=',q95_pos_7)
##
negPDF_7=NegPDF(7)
pyplot.hist(negPDF_7)
q05_neg_7=quantile(negPDF_7,0.05)
q95_neg_7=quantile(negPDF_7,0.95)
print('5% quantile of negPDF_7=',q05_neg_7)
print('95% quantile of negPDF_7=',q95_neg_7)

fig, axes = pyplot.subplots(1, 2)
fig.suptitle("Positive (left) and negative (right) PDF(7)")
axes[0].hist(posPDF_7)
axes[1].hist(negPDF_7)

##
posPDF_8=PosPDF(8)
pyplot.hist(posPDF_8)
q05_pos_8=quantile(posPDF_8,0.05)
q95_pos_8=quantile(posPDF_8,0.95)
print('5% quantile of posPDF_8=',q05_pos_8)
print('95% quantile of posPDF_8=',q95_pos_8)
##
negPDF_8=NegPDF(8)
pyplot.hist(negPDF_8)
q05_neg_8=quantile(negPDF_8,0.05)
q95_neg_8=quantile(negPDF_8,0.95)
print('5% quantile of negPDF_8=',q05_neg_8)
print('95% quantile of negPDF_8=',q95_neg_8)

fig, axes = pyplot.subplots(1, 2)
fig.suptitle("Positive (left) and negative (right) PDF(8)")
axes[0].hist(posPDF_8)
axes[1].hist(negPDF_8)

##########applications
##Agricultural land area_new
csv_data_agrland = pd.read_csv(fd.askopenfilename(),encoding = 'unicode_escape')
InputData_agrland = array(csv_data_agrland)
##Ammonia Emission
csv_data_AE = pd.read_csv(fd.askopenfilename(),encoding = 'unicode_escape')
InputData_AE = array(csv_data_AE)
##Energy use and biofuel production
csv_data_EUBP = pd.read_csv(fd.askopenfilename(),encoding = 'unicode_escape')
InputData_EUBP = array(csv_data_EUBP)
##Greenhouse gas emissions
csv_data_GGE = pd.read_csv(fd.askopenfilename(),encoding = 'unicode_escape')
InputData_GGE = array(csv_data_GGE)
## Farm Birds Index
csv_data_FBI = pd.read_csv(fd.askopenfilename(),encoding = 'unicode_escape')
InputData_FBI = array(csv_data_FBI)
####Tonnes/Hec
AE_over_agrland=InputData_AE/InputData_agrland
EUBP_over_agrland=InputData_EUBP/InputData_agrland
GGE_over_agrland=InputData_GGE/InputData_agrland
## 1st order difference
diff_AE_over_agrland=diff(AE_over_agrland)
diff_EUBP_over_agrland=diff(EUBP_over_agrland)
diff_GGE_over_agrland=diff(GGE_over_agrland)
diff_InputData_FBI=diff(InputData_FBI)
## country-series data (application2 in the paper)
num_country=shape(diff_AE_over_agrland)[0]
countrySeries=list([0 for i in range(num_country)])
coposCausMat_list=list([0 for i in range(num_country)])
conegCausMat_list=list([0 for i in range(num_country)])
index_max_posCaus=list([0 for i in range(num_country)])
index_min_negCaus=list([0 for i in range(num_country)])
for i in range(num_country):
  dif_ae_agrland=diff_AE_over_agrland[i]
  dif_eubp_agrland=diff_EUBP_over_agrland[i]
  dif_gge_agrland=diff_GGE_over_agrland[i]
  diff_id_fbi=diff_InputData_FBI[i]
  time_by_var_i=column_stack((dif_ae_agrland,dif_eubp_agrland
      ,dif_gge_agrland,diff_id_fbi))
  countrySeries[i]=time_by_var_i
  coposCausMat_list[i]=CoposCausMat(time_by_var_i)
  conegCausMat_list[i]=ConegCausMat(time_by_var_i)
  Altman1=coposCausMat_list[i]
  Altman2=coposCausMat_list[i]
  ##Altman_minus2=fill_diagonal(Altman,-2)
  Altman_dig_minus2 =Altman1;Altman_dig_2=Altman2
  Altman_dig_minus2[diag_indices_from(Altman_dig_minus2)]=-2
  index_max_posCaus[i]=unravel_index(argmax(Altman_dig_minus2, axis=None), Altman_dig_minus2.shape)
  Altman_dig_2[diag_indices_from(Altman_dig_2)]=3
  index_min_negCaus[i]=unravel_index(argmin(Altman_dig_2, axis=None), Altman_dig_2.shape)
## to latex table
print(tabulate(diff_AE_over_agrland, tablefmt="latex", floatfmt=".3f"))
print(tabulate(diff_EUBP_over_agrland, tablefmt="latex", floatfmt=".3f"))
print(tabulate(diff_GGE_over_agrland, tablefmt="latex", floatfmt=".3f"))
print(tabulate(diff_InputData_FBI, tablefmt="latex", floatfmt=".3f"))