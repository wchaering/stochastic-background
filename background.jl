using NoiseModel.jl

Tobs = 2.0^21
Nobs = 2^20

N = 100000 

model_curr = LISANoiseModel(4.5e-42, 4.5e-42, 4.5e-42, 4.5e-42, 9.2e-50, 9.2e-50, 9.2e-50, 9.2e-50, log(1e-11))
model_prop = LISANoiseModel(4.5e-42, 4.5e-42, 4.5e-42, 4.5e-42, 9.2e-50, 9.2e-50, 9.2e-50, 9.2e-50, log(1e-11))

XX0  =  0.3;
XX2  = -0.1005952381;
XX4  =  0.0140542328;
XX6  = -0.0010338567;
XX8  =  0.00004679056155;
XX10 = -0.000001434662976;
XXp=Array(Float64,N)
XXa=Array(Float64,N)
XX=Array(Float64,N)
f=Array(Float64,N)        #frequency
cosf=Array(Float64,N)     #for cosine terms in noise functions
sinf=Array(Float64,N)     #for sine terms in noise functions
SaScale=Array(Float64,N)  #acceleration noise component
Jumpscale = 0.05
LowXX = Array(Float64,N)
XXgw = Array(Float64,N)
Sgw = Array(Float64,N)
fstar = 0.00954269032
hconst = 2.26854354e-18

function CalcSgw(logOmega,f)
  ((3*hconst^2)/(4*pi^2*f^3))*exp(logOmega)
end

function LISANoise_XXp(sinf, model::LISANoiseModel)
  return 4.0*sinf*sinf*(model.Sp12 + model.Sp13 + model.Sp31 + model.Sp21)
end

function LISANoise_XXa(sinf, cosf, model::LISANoiseModel)
  return 16.0*sinf*sinf*(model.Sa12 + model.Sa13 + (model.Sa31 + model.Sa21)*cosf*cosf)
end

function LF_eqArmMich(fonfs)
  return XX0 + XX2*fonfs*fonfs + XX4*fonfs*fonfs*fonfs*fonfs + XX6*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs
  + XX8*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs + XX10*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs*fonfs;
end

for i=1:N
  f[i]=i/Tobs
  sinf[i] = sin(f[i]/fstar)
  cosf[i] = cos(f[i]/fstar)
  SaScale[i]=(2.0*pi*f[i])^(-4.0)*(1.0+(1.0e4*f[i])^(-2.0))
  XXp[i] = LISANoise_XXp(sinf[i], model_curr)
  XXa[i] = SaScale[i]*LISANoise_XXa(sinf[i], cosf[i], model_curr)
  Sgw[i] = CalcSgw(model_curr.LogO, f[i])
  LowXX[i] = LF_eqArmMich(f[i]/fstar)
  XXgw[i] = 4 * sinf[i]^2 * LowXX[i] * Sgw[i]
  XX[i] = (XXp[i]+XXa[i]) + XXgw[i]
end

noise_data = readdlm("CombinedData.dat")
data = noise_data[1:N]

logL_curr = LISANoise_Likelihood(data, XX, N)  #Calculate the log likelihood

COUNT = 100000
acc = 0.0

outfile=open("chains_combined.dat","w")
for m=1:COUNT
  if(m%100 == 0) println(m, ' ', acc/m) end
  model_prop = LISANoise_Proposals(model_curr, N)

  for i=1:N
    XXp[i] = LISANoise_XXp(sinf[i], model_prop)
    XXa[i] = SaScale[i]*LISANoise_XXa(sinf[i], cosf[i], model_prop)
    Sgw[i] = CalcSgw(model_prop.LogO,f[i])
    XXgw[i] = 4 * sinf[i]^2 * LowXX[i] * Sgw[i]
    XX[i] = (XXp[i] + XXa[i]) + XXgw[i]
  end

  logL_prop = LISANoise_Likelihood(data, XX, N)
  dlog = logL_prop - logL_curr
  hastings = exp(dlog)

  if hastings > rand()
    model_curr = model_prop
    logL_curr = logL_prop
    acc += 1.0
  end

  print(outfile, ' ', m, ' ',logL_curr, ' ', exp(model_curr.LogO), ' ', model_curr.Sp12/4.0e-42 - 1.0, ' ', (model_curr.Sp12+model_curr.Sp21+model_curr.Sp13+model_curr.Sp31)/4.0e-42-4.0, ' ',
  (model_curr.Sa12+model_curr.Sa21+model_curr.Sa13+model_curr.Sa31)/9.0e-50-4.0, "\n")

end
close(outfile)

plots = readdlm("chains_combined.dat")