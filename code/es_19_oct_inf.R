esercizio_informatici_19_oct = function(samples, min_unif = 0, max_unif = 1){

  beta_mixture = function(x){
    value = 0.3*(x**4) * (1-x) * (factorial(6)/(factorial(4))) + 0.7 * (x * (1-x)**7) * factorial(9) / factorial(7)
    return(value)
  }
  
  #x_i
  y_unif= runif(samples, min_unif, max_unif)
  print(y_unif)
  
  #f(x_i)
  f_x = c()
  print(f_x)
  #g(x_i)
  g_x = c()
  for(i in y_unif){
    f_x = c(f_x, beta_mixture(i))
    g_x = c(g_x, dunif(i, min_unif, max_unif))
  }
  print(f_x)
  print(g_x)
  #w_i = f(x_i)/g(x_i)
  w = f_x / g_x
  expected_value = sum(y_unif*w)/samples
  
  
  in_interval = c()
  for(f_x_i in f_x){
    in_interval = c(in_interval, f_x_i >= 0.45 && 0.55 >= f_x_i)
  }
  print(in_interval)
  pr = sum(in_interval) / length(in_interval)
  return(c(expected_value, pr))
}

