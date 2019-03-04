function s = AtimesB(r,z)
s = sum(sum(real(r).*real(z))) +  sum(sum(imag(r).*imag(z)));