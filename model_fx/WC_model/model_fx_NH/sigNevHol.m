function Fin = sigNevHol(in,Mn,Sn,Bn)
Fin = Mn./(1+(exp((-Sn.*in)./Mn)*(Mn-Bn))./Bn);