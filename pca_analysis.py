import numpy, csv, math
from ROOT import *
from numpy import linalg as LA
from project_data import project_data

plot_orig_dimensions 		= False
project_data_to_sig		= True
min_sig_eig_ratio		= -1
plot_projected_data_to_sig	= True
save_sig_eigenvectors_matrix	= True
project_data_to_bkg		= True
min_bkg_eig_ratio		= -1
plot_projected_data_to_bkg	= True
save_bkg_eigenvectors_matrix	= True
generate_test_to_sig		= True
generate_test_to_bkg		= True
test_file			= "../Higgs13TeV_test_118_130_ggh.csv"


reader = numpy.loadtxt("../Higgs13TeV_train_118_130_ggh.csv", delimiter=",")
X_train = reader[:,0:21].astype(float)
Y_train = reader[:,21]
nevents_train = len(Y_train)

#preparing data
ns = 0
nb = 0
for i in range(nevents_train):
  if(Y_train[i] == 1):
    ns += 1
  else:
    nb += 1

limit = ns
if(ns > nb) limit = nb
sX_train = [[]]
bX_train = [[]]
ns = 0
nb = 0
for i in range(nevents_train):
  if(Y_train[i] == 1 and ns < limit):
    sX_train.append( X_train[i] )
    ns += 1
  if(Y_train[i] == 0 and nb < limit):
    bX_train.append( X_train[i] )
    nb += 1

  

dimensions = len(X_train[0])
n_covariances = math.factorial(dimensions)/(math.factorial(dimensions-2)*2)
print "You have %i dimensions resulting to %i covariance matrices..." % (dimensions,n_covariances)


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
#the matrix to store signal covariances
sconv_matrix = [[0 for i in range(dimensions)] for j in range(dimensions)]
sconv_hist = TH2D("sconv_hist","",21,0,21,21,0,21)
#the matrix to store background covariances
bconv_matrix = [[0 for i in range(dimensions)] for j in range(dimensions)]
bconv_hist = TH2D("bconv_hist","",21,0,21,21,0,21)

#gets the sum of entries in each column
#vals_sum = numpy.sum(X_train, axis=0)
svals_sum = numpy.sum(sX_train, axis=0)
bvals_sum = numpy.sum(bX_train, axis=0)

#makes the PCA decomposition
for st_d in range(dimensions):
  print "Analysing through dimension: %i ..." % (st_d+1)
  for nd_d in range(dimensions):
    fcovariance = 0
    scovariance = 0
    bcovariance = 0
    #d1_mean = vals_sum[st_d]/float(nevents_train-1)
    d1_mean = svals_sum[st_d]/float(limit-1)
    for ip in range(limit):
      #d2_mean = vals_sum[nd_d]/float(nevents_train-1)
      d2_mean = bvals_sum[nd_d]/float(limit-1)
      fcovariance += (sX_train[ip][st_d]-d1_mean)*(bX_train[ip][nd_d]-d2_mean)
      #if(Y_train[ip] == 1):
	#scovariance += (X_train[ip][st_d]-d1_mean)*(X_train[ip][nd_d]-d2_mean)
      #else:
	#bcovariance += (X_train[ip][st_d]-d1_mean)*(X_train[ip][nd_d]-d2_mean)
    #full covariance matrix
    fconv_matrix[st_d][nd_d] = fcovariance/float(limit-1)
    fconv_hist.SetBinContent(nd_d+1,dimensions-st_d,fconv_matrix[st_d][nd_d])
    #signal matrix covariance
    #sconv_matrix[st_d][nd_d] = scovariance/float(nevents_train-1)
    #sconv_hist.SetBinContent(nd_d+1,dimensions-st_d,sconv_matrix[st_d][nd_d])
    #background covariance matrix
    #bconv_matrix[st_d][nd_d] = bcovariance/float(nevents_train-1)
    #bconv_hist.SetBinContent(nd_d+1,dimensions-st_d,bconv_matrix[st_d][nd_d])
  fconv_hist.GetXaxis().SetBinLabel(st_d+1,dim_names[st_d])
  fconv_hist.GetYaxis().SetBinLabel(dimensions-st_d,dim_names[st_d])
  #sconv_hist.GetXaxis().SetBinLabel(st_d+1,dim_names[st_d])
  #sconv_hist.GetYaxis().SetBinLabel(dimensions-st_d,dim_names[st_d])
  #bconv_hist.GetXaxis().SetBinLabel(st_d+1,dim_names[st_d])
  #bconv_hist.GetYaxis().SetBinLabel(dimensions-st_d,dim_names[st_d])
    

#exporting the covariance matrix in a txt file
full_conv_matrix = open('full_covariance_matrix.txt','w')
#sig_conv_matrix = open('sig_covariance_matrix.txt','w')
#bkg_conv_matrix = open('bkg_covariance_matrix.txt','w')
for i in range(dimensions):
  for j in range(dimensions):
    text1 = "%f" % fconv_matrix[i][j]
    full_conv_matrix.write( text1 )
    #text2 = "%f" % sconv_matrix[i][j]
    #sig_conv_matrix.write( text2 )
    #text3 = "%f" % bconv_matrix[i][j]
    #bkg_conv_matrix.write( text3 )
    if(j < dimensions-1):
      full_conv_matrix.write( "," )
      #sig_conv_matrix.write( "," )
      #bkg_conv_matrix.write( "," )
  full_conv_matrix.write( "\n" )
  #sig_conv_matrix.write( "\n" )
  #bkg_conv_matrix.write( "\n" )
#closing files
full_conv_matrix.close()
#sig_conv_matrix.close()
#bkg_conv_matrix.close()


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
#plot sig covariance
#sconv_hist.Draw("text")
#sconv_hist.SetTitle("Signal Covariance Matrix")
#gPad.Update()
#cv.Print("plots/Sig_covariance_matrix.png")
#plot bkg covariance
#bconv_hist.Draw("text")
#bconv_hist.SetTitle("Background Covariance Matrix")
#gPad.Update()
#cv.Print("plots/Bkg_covariance_matrix.png")


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

#find the signal eigenvectors and eigenvalues
#seigenvalues, seigenvectors = LA.eig(fconv_matrix)
#seigenvalues, seigenvectors = LA.eig(sconv_matrix)
#seig_vec_hist = TH2D("seig_vec_hist","",21,0,21,22,0,22)
#for i in range(len(seigenvalues)):
#  for j in range(len(seigenvectors[i])):
#    seig_vec_hist.SetBinContent(j+1,len(seigenvalues)-i,seigenvectors[i][j])
#  seig_vec_hist.SetBinContent(i+1,len(seigenvalues),seigenvalues[i])
#  seig_vec_hist.GetXaxis().SetBinLabel(i+1,dim_names[i])
#  seig_vec_hist.GetYaxis().SetBinLabel(len(seigenvalues)-i,dim_names[i])
#seig_vec_hist.GetYaxis().SetBinLabel(len(seigenvalues),"eigs")
#seig_vec_hist.Draw("text")
#seig_vec_hist.SetTitle("Signal Eigenvectors")
#cv.Print("plots/Sig_eigenvectors.png")

#find the background eigenvectors and eigenvalues
#beigenvalues, beigenvectors = LA.eig(bconv_matrix)
#beig_vec_hist = TH2D("beig_vec_hist","",21,0,21,21,0,21)
#for i in range(len(beigenvalues)):
#  for j in range(len(beigenvectors[i])):
#    beig_vec_hist.SetBinContent(j+1,len(beigenvalues)-i,beigenvectors[i][j])
#  beig_vec_hist.SetBinContent(i+1,len(beigenvalues),beigenvalues[i])
#  beig_vec_hist.GetXaxis().SetBinLabel(i+1,dim_names[i])
#  beig_vec_hist.GetYaxis().SetBinLabel(len(beigenvalues)-i,dim_names[i])
#beig_vec_hist.GetYaxis().SetBinLabel(len(beigenvalues),"eigs")
#beig_vec_hist.Draw("text")
#beig_vec_hist.SetTitle("Background Eigenvectors")
#cv.Print("plots/Bkg_eigenvectors.png")


#plot all pair dimensions and 1st, 2nd PCAs for signal and background
if(plot_orig_dimensions):
  cv = TCanvas("cv","",10,10,1000,1000)
  for st_d in range(dimensions):
    print "\nPlotting through dimension: %i ..." % (st_d+1)
    for nd_d in range(dimensions):
      gcv_s = TGraph()
      #gcv_b = TGraph()
      ipoint_s = 0
      #ipoint_b = 0
      d1_mean = svals_sum[st_d]/float(limit-1)
      for ip in range(limit):
	d2_mean = bvals_sum[nd_d]/float(limit-1)
	#if(Y_train[ip] == 1):
	gcv_s.SetPoint(ipoint_s,sX_train[ip][st_d]-d1_mean,bX_train[ip][nd_d]-d2_mean)
	ipoint_s += 1
	#else:
	  #gcv_b.SetPoint(ipoint_b,X_train[ip][st_d]-d1_mean,X_train[ip][nd_d]-d2_mean)
	  #ipoint_b += 1
      gcv_s.SetMarkerColor(kBlue)
      gcv_s.SetMarkerStyle(2) #6
      #gcv_b.SetMarkerColor(kRed)
      #gcv_b.SetMarkerStyle(2)
      mg = TMultiGraph()
      mg.Add(gcv_s)
      #mg.Add(gcv_b)
      mg.Draw("ap")
      mg.GetXaxis().SetTitle(dim_names[st_d])
      mg.GetYaxis().SetTitle(dim_names[nd_d])
      leg = TLegend(0.75,0.7,0.95,0.92)
      leg.SetHeader("Original Space")
      #leg.AddEntry(gcv_s,"Signal","p")
      #leg.AddEntry(gcv_b,"Background","p")
      f_eig_1 = feigenvectors[st_d][st_d]/feigenvectors[nd_d][st_d]
      f_eig_2 = feigenvectors[st_d][nd_d]/feigenvectors[nd_d][nd_d]
      #s_eig_1 = seigenvectors[st_d][st_d]/seigenvectors[nd_d][st_d]
      #s_eig_2 = seigenvectors[st_d][nd_d]/seigenvectors[nd_d][nd_d]
      #b_eig_1 = beigenvectors[st_d][st_d]/beigenvectors[nd_d][st_d]
      #b_eig_2 = beigenvectors[st_d][nd_d]/beigenvectors[nd_d][nd_d]
      #sorting the PCAs
      feig1f_par = 0
      feig2f_par = 0
      #seig1f_par = 0
      #seig2f_par = 0
      #beig1f_par = 0
      #beig2f_par = 0
      if(feigenvalues[st_d] > feigenvalues[nd_d]):
	feig1f_par = f_eig_1
	feig2f_par = f_eig_2
      else:
	feig1f_par = f_eig_2
	feig2f_par = f_eig_1
      #if(seigenvalues[st_d] > seigenvalues[nd_d]):
	#seig1f_par = s_eig_1
	#seig2f_par = s_eig_2
      #else:
	#seig1f_par = s_eig_2
	#seig2f_par = s_eig_1
      #if(beigenvalues[st_d] > beigenvalues[nd_d]):
	#beig1f_par = b_eig_1
	#beig2f_par = b_eig_2
      #else:
	#beig1f_par = b_eig_2
	#beig2f_par = b_eig_1
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
      #seig1f = TF1("seig1f","[0]*x")
      #seig1f.SetParameter(0,seig1f_par)
      #lim = math.fabs(seig1f.GetX(fy_max,-1.e3,1.e3))
      #seig1f.SetRange(-lim,lim)
      #seig1f.SetLineColor(kOrange+2)
      #seig1f.SetLineWidth(8)
      #seig2f = TF1("seig2f","[0]*x")
      #seig2f.SetParameter(0,seig2f_par)
      #lim = math.fabs(seig2f.GetX(fy_max,-1.e3,1.e3))
      #seig2f.SetRange(-lim,lim)
      #seig2f.SetLineColor(kOrange-3)
      #seig2f.SetLineWidth(6)
      #beig1f = TF1("beig1f","[0]*x")
      #beig1f.SetParameter(0,beig1f_par)
      #lim = math.fabs(beig1f.GetX(fy_max,-1.e3,1.e3))
      #beig1f.SetRange(-lim,lim)
      #beig1f.SetLineColor(kYellow+2)
      #beig1f.SetLineWidth(4)
      #beig2f = TF1("beig2f","[0]*x")
      #beig2f.SetParameter(0,beig2f_par)
      #lim = math.fabs(beig2f.GetX(fy_max,-1.e3,1.e3))
      #beig2f.SetRange(-lim,lim)
      #beig2f.SetLineColor(kYellow-3)
      #beig2f.SetLineWidth(2)
      feig1f.Draw("l,same")
      feig2f.Draw("l,same")
      #seig1f.Draw("l,same")
      #seig2f.Draw("l,same")
      #beig1f.Draw("l,same")
      #beig2f.Draw("l,same")
      leg.AddEntry(feig1f,"1^{st} PCA","l")
      leg.AddEntry(feig2f,"2^{nd} PCA","l")
      #leg.AddEntry(seig1f,"Sig 1^{st} PCA","l")
      #leg.AddEntry(seig2f,"Sig 2^{nd} PCA","l")
      #leg.AddEntry(beig1f,"Bkg 1^{st} PCA","l")
      #leg.AddEntry(beig2f,"Bkg 2^{nd} PCA","l")
      #leg.AddEntry(seig1f,("%fx" % seig1f_par),"l")
      #leg.AddEntry(seig2f,("%fx" % seig2f_par),"l")
      #leg.AddEntry(beig1f,("%fx" % beig1f_par),"l")
      #leg.AddEntry(beig2f,("%fx" % beig2f_par),"l")
      leg.Draw()
      gPad.Update()
      name = "plots/Covariance_d%i_d%i.png" % (st_d,nd_d)
      cv.Print(name)

cv.Close()
cv = TCanvas("cv","",10,10,1000,1000)

#transposing the original dataset
if(project_data_to_sig or project_data_to_bkg):
  print "\nTransposing original dataset..."
  transp_dataset = [[0 for i in range(nevents_train)] for j in range(dimensions)]
  for i in range(dimensions):
    for j in range(nevents_train):
      transp_dataset[i][j] = X_train[j][i]

  #transposing the signal eigenvectors
  if(project_data_to_sig):
    print "Transposing signal eigenvectors..."
    transp_sig_eigenvectors = [[0 for i in range(dimensions)] for j in range(dimensions)]
    for i in range(dimensions):
      for j in range(dimensions):
	transp_sig_eigenvectors[i][j] = seigenvectors[j][i]

    #sort the eigenvectors based on the eigenvalues (highest to smallest - principal to sub-principal components)
    print "Sorting signal eigenvectors..."
    sorted_indexes = numpy.argsort(seigenvalues)
    decreasing = list(reversed(sorted_indexes))
    sorted_sig_eigenvectors = []
    for i in range(dimensions):
      index = decreasing[i]
      #remove sub-principal components
      if(seigenvalues[index]/numpy.amax(seigenvalues) < min_sig_eig_ratio and min_sig_eig_ratio != -1):
	continue
      sorted_sig_eigenvectors.append( transp_sig_eigenvectors[index] )
    transp_sig_eigenvectors = sorted_sig_eigenvectors
    
    #save signal eigenvectors
    if(save_sig_eigenvectors_matrix):
      sig_eigenvectors_file = open('signal_eigenvectors.txt','w')
      for i in range(len(transp_sig_eigenvectors)):
	for j in range(len(transp_sig_eigenvectors[0])):
	  info = "%.15f" % transp_sig_eigenvectors[i][j]
	  sig_eigenvectors_file.write( info )
	  if(j < len(transp_sig_eigenvectors[0])-1):
	    sig_eigenvectors_file.write( "," )
	sig_eigenvectors_file.write( "\n" )
      sig_eigenvectors_file.close()
      

    #projecting the original dataset into the signal eigenvectors phase space
    print "Projecting dataset through signal eigenvectors space..."
    projected_dataset_to_sig = [[0 for j in range(len(transp_sig_eigenvectors))] for j in range(nevents_train)]
    for irow in range(len(transp_sig_eigenvectors)):
      for ientry in range(nevents_train):
	entry = 0
	for icol in range(len(transp_sig_eigenvectors[0])):
	  #FinalData = RowFeatureVector x RowDataAdjusted					#data mean
	  entry += transp_sig_eigenvectors[irow][icol]*(transp_dataset[icol][ientry] - vals_sum[icol]/float(nevents_train-1))
	projected_dataset_to_sig[ientry][irow] = entry

    #plotting the projected data
    if(plot_projected_data_to_sig):
      for st_d in range(len(projected_dataset_to_sig[0])):
	print "\nPlotting through dimension: %i ..." % (st_d+1)
	for nd_d in range(len(projected_dataset_to_sig[0])):
	  gcv_s2 = TGraph()
	  gcv_b2 = TGraph()
	  ipoint_s = 0
	  ipoint_b = 0
	  for ip in range(nevents_train):
	    if(Y_train[ip] == 1):
	      gcv_s2.SetPoint(ipoint_s,projected_dataset_to_sig[ip][st_d],projected_dataset_to_sig[ip][nd_d])
	      ipoint_s += 1
	    else:
	      gcv_b2.SetPoint(ipoint_b,projected_dataset_to_sig[ip][st_d],projected_dataset_to_sig[ip][nd_d])
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
	  name = "projected_plots_sig_space/Covariance_d%i_d%i.png" % (st_d,nd_d)
	  cv.Print(name)

    if(generate_test_to_sig):
      to_test = numpy.loadtxt(test_file, delimiter=",")
      X_test = reader[:,0:21].astype(float)
      nevents_test = len(X_test)
      pX_test = project_data(X_test,'signal_eigenvectors.txt')

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
	  name = "test_on_sig_space/Covariance_d%i_d%i.png" % (st_d,nd_d)
	  cv.Print(name)



  #transposing the background eigenvectors
  if(project_data_to_bkg):
    print "\nTransposing background eigenvectors..."
    transp_bkg_eigenvectors = [[0 for i in range(dimensions)] for j in range(dimensions)]
    for i in range(dimensions):
      for j in range(dimensions):
	transp_bkg_eigenvectors[i][j] = beigenvectors[j][i]

    #sort the eigenvectors based on the eigenvalues (highest to smallest - principal to sub-principal components)
    print "Sorting background eigenvectors..."
    sorted_indexes = numpy.argsort(beigenvalues)
    decreasing = list(reversed(sorted_indexes))
    sorted_bkg_eigenvectors = []
    for i in range(dimensions):
      index = decreasing[i]
      #remove sub-principal components
      if(beigenvalues[index]/numpy.amax(beigenvalues) < min_bkg_eig_ratio and min_bkg_eig_ratio != -1):
	continue
      sorted_bkg_eigenvectors.append( transp_bkg_eigenvectors[index] )
    transp_bkg_eigenvectors = sorted_bkg_eigenvectors
    
    #save background eigenvectors
    if(save_bkg_eigenvectors_matrix):
      bkg_eigenvectors_file = open('background_eigenvectors.txt','w')
      for i in range(len(transp_bkg_eigenvectors)):
	for j in range(len(transp_bkg_eigenvectors[0])):
	  info = "%.15f" % transp_bkg_eigenvectors[i][j]
	  bkg_eigenvectors_file.write( info )
	  if(j < len(transp_bkg_eigenvectors[0])-1):
	    bkg_eigenvectors_file.write( "," )
	bkg_eigenvectors_file.write( "\n" )
      bkg_eigenvectors_file.close()      
      

    #projecting the original dataset into the signal eigenvectors phase space
    print "Projecting dataset through background eigenvectors space..."
    projected_dataset_to_bkg = [[0 for j in range(len(transp_bkg_eigenvectors))] for j in range(nevents_train)]
    for irow in range(len(transp_bkg_eigenvectors)):
      for ientry in range(nevents_train):
	entry = 0
	for icol in range(dimensions):
	  #FinalData = RowFeatureVector x RowDataAdjusted					#data mean
	  entry += transp_bkg_eigenvectors[irow][icol]*(transp_dataset[icol][ientry] - vals_sum[icol]/float(nevents_train-1))
	projected_dataset_to_bkg[ientry][irow] = entry

    #plotting the projected data
    if(plot_projected_data_to_bkg):
      for st_d in range(len(projected_dataset_to_bkg[0])):
	print "\nPlotting through dimension: %i ..." % (st_d+1)
	for nd_d in range(len(projected_dataset_to_bkg[0])):
	  gcv_s3 = TGraph()
	  gcv_b3 = TGraph()
	  ipoint_s = 0
	  ipoint_b = 0
	  for ip in range(nevents_train):
	    if(Y_train[ip] == 1):
	      gcv_s3.SetPoint(ipoint_s,projected_dataset_to_bkg[ip][st_d],projected_dataset_to_bkg[ip][nd_d])
	      ipoint_s += 1
	    else:
	      gcv_b3.SetPoint(ipoint_b,projected_dataset_to_bkg[ip][st_d],projected_dataset_to_bkg[ip][nd_d])
	      ipoint_b += 1
	  gcv_s3.SetMarkerColor(kBlue)
	  gcv_s3.SetMarkerStyle(2) #6
	  gcv_b3.SetMarkerColor(kRed)
	  gcv_b3.SetMarkerStyle(2)
	  mg = TMultiGraph()
	  mg.Add(gcv_s3)
	  mg.Add(gcv_b3)
	  mg.Draw("ap")
	  mg.GetXaxis().SetTitle("#bf{#hat{e}}_{%i}" % st_d)
	  mg.GetYaxis().SetTitle("#bf{#hat{e}}_{%i}" % nd_d)
	  leg = TLegend(0.15,0.8,0.37,0.92)
	  leg.SetHeader("Bkg Space")
	  leg.AddEntry(gcv_s3,"Signal","p")
	  leg.AddEntry(gcv_b3,"Background","p")
	  leg.Draw()
	  gPad.Update()
	  name = "projected_plots_bkg_space/Covariance_d%i_d%i.png" % (st_d,nd_d)
	  cv.Print(name)
	  
	  
    if(generate_test_to_bkg):
      to_test = numpy.loadtxt(test_file, delimiter=",")
      X_test = reader[:,0:21].astype(float)
      nevents_test = len(X_test)
      pX_test = project_data(X_test,'background_eigenvectors.txt')

      for st_d in range(len(pX_test[0])):
	print "\nPlotting through dimension: %i ..." % (st_d+1)
	for nd_d in range(len(pX_test[0])):
	  gcv_s5 = TGraph()
	  gcv_b5 = TGraph()
	  ipoint_s = 0
	  ipoint_b = 0
	  for ip in range(nevents_train):
	    if(Y_train[ip] == 1):
	      gcv_s5.SetPoint(ipoint_s,pX_test[ip][st_d],pX_test[ip][nd_d])
	      ipoint_s += 1
	    else:
	      gcv_b5.SetPoint(ipoint_b,pX_test[ip][st_d],pX_test[ip][nd_d])
	      ipoint_b += 1
	  gcv_s5.SetMarkerColor(kBlue)
	  gcv_s5.SetMarkerStyle(2) #6
	  gcv_b5.SetMarkerColor(kRed)
	  gcv_b5.SetMarkerStyle(2)
	  mg = TMultiGraph()
	  mg.Add(gcv_s5)
	  mg.Add(gcv_b5)
	  mg.Draw("ap")
	  mg.GetXaxis().SetTitle("#bf{#hat{e}}_{%i}" % st_d)
	  mg.GetYaxis().SetTitle("#bf{#hat{e}}_{%i}" % nd_d)
	  leg = TLegend(0.15,0.8,0.37,0.92)
	  leg.SetHeader("Bkg Space")
	  leg.AddEntry(gcv_s5,"Signal","p")
	  leg.AddEntry(gcv_b5,"Background","p")
	  leg.Draw()
	  gPad.Update()
	  name = "test_on_bkg_space/Covariance_d%i_d%i.png" % (st_d,nd_d)
	  cv.Print(name)
      