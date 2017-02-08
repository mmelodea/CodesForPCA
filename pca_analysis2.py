import numpy, csv, math
from ROOT import *
from numpy import linalg as LA
from project_data import project_data

test_file			= "../Higgs13TeV_test_118_130_ggh.csv"

reader = numpy.loadtxt("../Higgs13TeV_train_118_130_ggh.csv", delimiter=",")
X_train = reader[:,0:21].astype(float)
Y_train = reader[:,21]
nevents_train = len(Y_train)
dimensions = len(X_train[0])
n_covariances = math.factorial(dimensions)/(math.factorial(dimensions-2)*2)
print "You have %i dimensions resulting to %i covariance matrices..." % (dimensions,n_covariances)

#preparing data
ns = 0
nb = 0
for i in range(nevents_train):
  if(Y_train[i] == 1):
    ns += 1
  else:
    nb += 1

limit = ns
if(ns > nb):
  limit = nb
sX_train = [[0 for i in range(dimensions)] for j in range(limit)]
bX_train = [[0 for i in range(dimensions)] for j in range(limit)]
ns = 0
nb = 0
for i in range(nevents_train):
  if(Y_train[i] == 1 and ns < limit):
    for j in range(dimensions):
      sX_train[ns][j] = X_train[ns][j]
    ns += 1
  if(Y_train[i] == 0 and nb < limit):
    for j in range(dimensions):
      bX_train[nb][j] = X_train[nb][j]
    nb += 1
  

dim_names = ["l1 p_{T}","l1 #eta","l1 #phi",
	     "l2 p_{T}","l2 #eta","l2 #phi",
	     "l3 p_{T}","l3 #eta","l3 #phi",
	     "l4 p_{T}","l4 #eta","l4 #phi",
	     "j1 p_{T}","j1 #eta","j1 #phi","j1 E",
	     "j2 p_{T}","j2 #eta","j2 #phi","j2 E",
	     "Njets"]

#the matrix to store all covariances
fconv_matrix = [[0 for i in range(dimensions)] for j in range(dimensions)]
fconv_hist = TH2D("fconv_hist","",21,0,21,21,0,21)

#gets the sum of entries in each column
svals_sum = numpy.sum(sX_train, axis=0)
bvals_sum = numpy.sum(bX_train, axis=0)

#makes the PCA decomposition
for st_d in range(dimensions):
  print "Analysing through dimension: %i ..." % (st_d+1)
  for nd_d in range(dimensions):
    fcovariance = 0
    d1_mean = svals_sum[st_d]/float(limit-1)
    for ip in range(limit):
      d2_mean = bvals_sum[nd_d]/float(limit-1)
      fcovariance += (sX_train[ip][st_d]-d1_mean)*(bX_train[ip][nd_d]-d2_mean)
    fconv_matrix[st_d][nd_d] = fcovariance/float(limit-1)
    fconv_hist.SetBinContent(nd_d+1,dimensions-st_d,fconv_matrix[st_d][nd_d])
  fconv_hist.GetXaxis().SetBinLabel(st_d+1,dim_names[st_d])
  fconv_hist.GetYaxis().SetBinLabel(dimensions-st_d,dim_names[st_d])    

#exporting the covariance matrix in a txt file
full_conv_matrix = open('full_covariance_matrix.txt','w')
for i in range(dimensions):
  for j in range(dimensions):
    text1 = "%f" % fconv_matrix[i][j]
    full_conv_matrix.write( text1 )
    if(j < dimensions-1):
      full_conv_matrix.write( "," )
  full_conv_matrix.write( "\n" )
#closing files
full_conv_matrix.close()


gROOT.SetBatch()
gStyle.SetOptStat(0)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetTitleOffset(1.4,"y")
gStyle.SetPaintTextFormat("0.2e")
cv = TCanvas("cv","",10,10,2300,1000)
#plot full covariance
fconv_hist.Draw("text")
fconv_hist.SetTitle("Full Covariance Matrix")
gPad.Update()
cv.Print("plots/Full_covariance_matrix.png")

feigenvalues, feigenvectors = LA.eig(fconv_matrix)
feig_vec_hist = TH2D("feig_vec_hist","",21,0,21,22,0,22)
for i in range(len(feigenvalues)):
  for j in range(len(feigenvectors[i])):
    feig_vec_hist.SetBinContent(j+1,len(feigenvalues)-i,feigenvectors[i][j])
  feig_vec_hist.SetBinContent(i+1,len(feigenvalues),feigenvalues[i])
  feig_vec_hist.GetXaxis().SetBinLabel(i+1,dim_names[i])
  feig_vec_hist.GetYaxis().SetBinLabel(len(feigenvalues)-i,dim_names[i])
feig_vec_hist.GetYaxis().SetBinLabel(len(feigenvalues),"eigs")
feig_vec_hist.Draw("text")
feig_vec_hist.SetTitle("Full Eigenvectors")
cv.Print("plots/Full_eigenvectors.png")


#plot all pair dimensions and 1st, 2nd PCAs for signal and background
cv = TCanvas("cv","",10,10,1000,1000)
for st_d in range(dimensions):
  print "\nPlotting through dimension: %i ..." % (st_d+1)
  for nd_d in range(dimensions):
    gcv_s = TGraph()
    ipoint_s = 0
    d1_mean = svals_sum[st_d]/float(limit-1)
    for ip in range(limit):
      d2_mean = bvals_sum[nd_d]/float(limit-1)
      gcv_s.SetPoint(ipoint_s,sX_train[ip][st_d]-d1_mean,bX_train[ip][nd_d]-d2_mean)
      ipoint_s += 1
    gcv_s.SetMarkerColor(kBlue)
    gcv_s.SetMarkerStyle(2) #6
    mg = TMultiGraph()
    mg.Add(gcv_s)
    mg.Draw("ap")
    mg.GetXaxis().SetTitle(dim_names[st_d])
    mg.GetYaxis().SetTitle(dim_names[nd_d])
    leg = TLegend(0.75,0.7,0.95,0.92)
    leg.SetHeader("Original Space")
    f_eig_1 = feigenvectors[st_d][st_d]/feigenvectors[nd_d][st_d]
    f_eig_2 = feigenvectors[st_d][nd_d]/feigenvectors[nd_d][nd_d]
    #sorting the PCAs
    feig1f_par = 0
    feig2f_par = 0
    if(feigenvalues[st_d] > feigenvalues[nd_d]):
      feig1f_par = f_eig_1
      feig2f_par = f_eig_2
    else:
      feig1f_par = f_eig_2
      feig2f_par = f_eig_1
    fy_max = numpy.max(sX_train,axis=1)[nd_d]
    feig1f = TF1("feig1f","[0]*x")
    feig1f.SetParameter(0,feig1f_par)
    lim = math.fabs(feig1f.GetX(fy_max,-1.e3,1.e3))
    feig1f.SetRange(-lim,lim)
    feig1f.SetLineColor(kOrange+2)
    feig1f.SetLineWidth(8)
    feig2f = TF1("feig2f","[0]*x")
    feig2f.SetParameter(0,feig2f_par)
    lim = math.fabs(feig2f.GetX(fy_max,-1.e3,1.e3))
    feig1f.Draw("l,same")
    feig2f.Draw("l,same")
    leg.AddEntry(feig1f,"1^{st} PCA","l")
    leg.AddEntry(feig2f,"2^{nd} PCA","l")
    leg.Draw()
    gPad.Update()
    name = "plots/Covariance_d%i_d%i.png" % (st_d,nd_d)
    cv.Print(name)

cv.Close()
cv = TCanvas("cv","",10,10,1000,1000)

#transposing the original dataset
print "\nTransposing original dataset..."
transp_dataset = [[0 for i in range(nevents_train)] for j in range(dimensions)]
for i in range(dimensions):
  for j in range(nevents_train):
    transp_dataset[i][j] = X_train[j][i]

#transposing the signal eigenvectors
print "Transposing signal eigenvectors..."
transp_eigenvectors = [[0 for i in range(dimensions)] for j in range(dimensions)]
for i in range(dimensions):
  for j in range(dimensions):
    transp_eigenvectors[i][j] = feigenvectors[j][i]

#sort the eigenvectors based on the eigenvalues (highest to smallest - principal to sub-principal components)
print "Sorting signal eigenvectors..."
sorted_indexes = numpy.argsort(feigenvalues)
decreasing = list(reversed(sorted_indexes))
sorted_eigenvectors = []
for i in range(dimensions):
  index = decreasing[i]
  sorted_eigenvectors.append( transp_eigenvectors[index] )
transp_eigenvectors = sorted_eigenvectors
    
#save signal eigenvectors
eigenvectors_file = open('full_eigenvectors.txt','w')
for i in range(len(transp_eigenvectors)):
  for j in range(len(transp_eigenvectors[0])):
    info = "%.15f" % transp_eigenvectors[i][j]
    eigenvectors_file.write( info )
    if(j < len(transp_eigenvectors[0])-1):
      eigenvectors_file.write( "," )
  eigenvectors_file.write( "\n" )
eigenvectors_file.close()
      

#projecting the original dataset into the signal eigenvectors phase space
print "Projecting dataset through signal eigenvectors space..."
projected_dataset = [[0 for j in range(len(transp_eigenvectors))] for j in range(nevents_train)]
for irow in range(len(transp_eigenvectors)):
  for ientry in range(nevents_train):
    entry = 0
    for icol in range(len(transp_eigenvectors[0])):
      data_mean = 0
      if(Y_train[ientry] == 1):
	data_mean = svals_sum[icol]/float(limit-1)
      else:
	data_mean = bvals_sum[icol]/float(limit-1)
      #FinalData = RowFeatureVector x RowDataAdjusted
      entry += transp_eigenvectors[irow][icol]*(transp_dataset[icol][ientry] - data_mean)
    projected_dataset[ientry][irow] = entry

for st_d in range(len(projected_dataset[0])):
  print "\nPlotting through dimension: %i ..." % (st_d+1)
  for nd_d in range(len(projected_dataset[0])):
    gcv_s2 = TGraph()
    gcv_b2 = TGraph()
    ipoint_s = 0
    ipoint_b = 0
    for ip in range(nevents_train):
      if(Y_train[ip] == 1):
	gcv_s2.SetPoint(ipoint_s,projected_dataset[ip][st_d],projected_dataset[ip][nd_d])
	ipoint_s += 1
      else:
	gcv_b2.SetPoint(ipoint_b,projected_dataset[ip][st_d],projected_dataset[ip][nd_d])
	ipoint_b += 1
    gcv_s2.SetMarkerColor(kBlue)
    gcv_s2.SetMarkerStyle(2) #6
    gcv_b2.SetMarkerColor(kRed)
    gcv_b2.SetMarkerStyle(2)
    mg = TMultiGraph()
    mg.Add(gcv_s2)
    mg.Add(gcv_b2)
    mg.Draw("ap")
    mg.GetXaxis().SetTitle("#bf{#hat{e}}_{%i}" % st_d)
    mg.GetYaxis().SetTitle("#bf{#hat{e}}_{%i}" % nd_d)
    leg = TLegend(0.15,0.8,0.37,0.92)
    leg.SetHeader("Sig Space")
    leg.AddEntry(gcv_s2,"Signal","p")
    leg.AddEntry(gcv_b2,"Background","p")
    leg.Draw()
    gPad.Update()
    name = "projected_plots/Covariance_d%i_d%i.png" % (st_d,nd_d)
    cv.Print(name)


to_test = numpy.loadtxt(test_file, delimiter=",")
X_test = reader[:,0:21].astype(float)
nevents_test = len(X_test)
pX_test = project_data(X_test,'full_eigenvectors.txt')

for st_d in range(len(pX_test[0])):
  print "\nPlotting through dimension: %i ..." % (st_d+1)
  for nd_d in range(len(pX_test[0])):
    gcv_s4 = TGraph()
    gcv_b4 = TGraph()
    ipoint_s = 0
    ipoint_b = 0
    for ip in range(nevents_train):
      if(Y_train[ip] == 1):
	gcv_s4.SetPoint(ipoint_s,pX_test[ip][st_d],pX_test[ip][nd_d])
	ipoint_s += 1
      else:
	gcv_b4.SetPoint(ipoint_b,pX_test[ip][st_d],pX_test[ip][nd_d])
	ipoint_b += 1
    gcv_s4.SetMarkerColor(kBlue)
    gcv_s4.SetMarkerStyle(2) #6
    gcv_b4.SetMarkerColor(kRed)
    gcv_b4.SetMarkerStyle(2)
    mg = TMultiGraph()
    mg.Add(gcv_s4)
    mg.Add(gcv_b4)
    mg.Draw("ap")
    mg.GetXaxis().SetTitle("#bf{#hat{e}}_{%i}" % st_d)
    mg.GetYaxis().SetTitle("#bf{#hat{e}}_{%i}" % nd_d)
    leg = TLegend(0.15,0.8,0.37,0.92)
    leg.SetHeader("Bkg Space")
    leg.AddEntry(gcv_s4,"Signal","p")
    leg.AddEntry(gcv_b4,"Background","p")
    leg.Draw()
    gPad.Update()
    name = "new_test/Covariance_d%i_d%i.png" % (st_d,nd_d)
    cv.Print(name)      