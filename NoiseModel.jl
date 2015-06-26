type LISANoiseModel
  Sp12::Float64
  Sp21::Float64
  Sp13::Float64
  Sp31::Float64
  Sa12::Float64
  Sa21::Float64
  Sa13::Float64
  Sa31::Float64
  LogO::Float64
end

function LISANoise_XXp(sinf, model::LISANoiseModel)
  return 4.0*sinf*sinf*(model.Sp12 + model.Sp13 + model.Sp31 + model.Sp21)
end

function LISANoise_XXa(sinf, cosf, model::LISANoiseModel)
  return 16.0*sinf*sinf*(model.Sa12 + model.Sa13 + (model.Sa31 + model.Sa21)*cosf*cosf)
end

function LISANoise_Proposals(curr::LISANoiseModel, N)
  L=5.0e9
  Sshot = 4.0e-22/(4.0*L*L);  #Position noise
  Sa = 9.0e-30/(4.0*L*L);     #acceleration noise
  sqN = sqrt(N)/10.0;
  #sqD = sqrt(1.0/heat);
  p12=Float64
  m12=Float64
  prop = LISANoiseModel(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,log(1e-11))

  while ( (prop.Sp12/Sshot > 10.0) ||  (prop.Sp12/Sshot < 0.1)  || (prop.Sp21/Sshot > 10.0) ||  (prop.Sp21/Sshot < 0.1));
    p12 =  Sshot*randn()/sqN;
    m12 =  10.0*Sshot*randn()/sqN;
    prop.Sp12 = curr.Sp12 + (p12+m12)/2.0;
    prop.Sp21 = curr.Sp21 + (p12-m12)/2.0;
  end

    while ( (prop.Sp13/Sshot > 10.0) ||  (prop.Sp13/Sshot < 0.1)  || (prop.Sp31/Sshot > 10.0) ||  (prop.Sp31/Sshot < 0.1));
        p12 =  Sshot*randn()/sqN;
        m12 =  10.0*Sshot*randn()/sqN;
        prop.Sp13 = curr.Sp13 + (p12+m12)/2.0;
        prop.Sp31 = curr.Sp31 + (p12-m12)/2.0;
    end

    while ( (prop.Sa12/Sa > 10.0) ||  (prop.Sa12/Sa < 0.1)  || (prop.Sa21/Sa > 10.0) ||  (prop.Sa21/Sa < 0.1));
        p12 =  Sa*randn()/sqN;
        m12 =  10.0*Sa*randn()/sqN;
        prop.Sa12 = curr.Sa12 + (p12+m12)/2.0;
        prop.Sa21 = curr.Sa21 + (p12-m12)/2.0;
    end

    while ( (prop.Sa13/Sa > 10.0) ||  (prop.Sa13/Sa < 0.1)  || (prop.Sa31/Sa > 10.0) ||  (prop.Sa31/Sa < 0.1));
        p12 =  Sa*randn()/sqN;
      m12 =  10.0*Sa*randn()/sqN;
     prop.Sa13 = curr.Sa13 + (p12+m12)/2.0;
      prop.Sa31 = curr.Sa31 + (p12-m12)/2.0;
    end
    prop.LogO = curr.LogO + (randn()*Jumpscale)

  return prop;
end

function LISANoise_Likelihood(data, model, N)
  detx=0.0;
  coef = Array(Float64,N)
  for i=1:N
    detx += log(model[i]);
    coef[i] = 1.0/model[i];
  end

  product = 0.0;
  for i=1:N
    product += coef[i]*data[i];
  end
  return -0.5*(2.0*product + 2.0*detx);
end