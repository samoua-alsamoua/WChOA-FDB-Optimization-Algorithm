function [ fitness ] = testFunction( x, fhd, fNumber )
          
        fitness = feval(fhd, x, fNumber);
    
end
