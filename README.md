This is code for ASP Vision (CVPR 2016) paper. It is built out of the MatConvNetbeta16 distribution, hence the extra code. Please look in the Examples folder for the ASP scripts. 

Our additions in Examples:
- Gw.m - file for Gabor wavelet creation for ASPs
- Amp.m - create Gaussian envelope for ASPs (used by GW.m)
- loadASPparameters.m - creates the tile of ASPs similar to the hardware setup
- cnn train.m has been modified for our purposes

Since Matconvnet is built on top of GPU computing, you will need to compile Matconvnet before you can use it. Please refer to the website below for more details (especially for the deprecated beta16 version). 

Note: This is research quality code (so probably will have bugs/lacks of detailed comments). Please contact us for further questions. We currently do not have plans to update our code for the newest version of Matconvnet or other deep learning packages. 

Please see our website: http://molnargroup.ece.cornell.edu/people/suren-jayasuriya/asp-vision/ for more details. Access to real ASP data is only granted on a case-by-case basis for educational/research purposes only, please contact Al Molnar or Suren Jayasuriya for your request.



# MatConvNet: CNNs for MATLAB

**MatConvNet** is a MATLAB toolbox implementing *Convolutional Neural
Networks* (CNNs) for computer vision applications. It is simple,
efficient, and can run and learn state-of-the-art CNNs. Several
example CNNs are included to classify and encode images. Please visit
the [homepage](http://www.vlfeat.org/matconvnet) to know more.
