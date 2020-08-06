
#
# Evaluate the compliance
#
function Compliance(U,F)

    return dot(U,F)

end

#
# Sensitivity analysis - Compliance
#
function DCompliance(ne,Ke,coord,connect,U,ρ,simp)

    # Initialize the output array
    d = zeros(ne)

    # Valor mínimo do simp
    ρmin = 1E-3

    # Aloca antes do loop para dar uma acelerada
    ul = @MVector zeros(8)

    # Loop
    for m=1:ne

         # Obtem os gls do elemento
         dofs = Dofs_ele(m,connect)

         # Extract local displacements (4x1)
         ul .= U[dofs]

         # Sensitivity w.r.t Ae
         derivada_simp = ((simp*(1.0-ρmin))*(ρmin+ρ[m]*(1-ρmin))^(simp-1.0))
         d[m] = -derivada_simp*dot(ul,Ke*ul)

    end #ele

   return d

end

#
# Volume. Here, volume stands for the sum of all the design variables
#
function Volume(ne::Int64,ρ::Array{Float64})

     sum(ρ)/ne

end

#
# Sensivity analysis - Volume
#
function dVolume(ne::Int64)
        ones(ne)/ne
end


 #
 # Optimality Criteria (VL the volume limit)
 #
 function CO_Flex(ne,x,dFlex,dVol,coord,connect,VL)


     # Basic setup for the algorithm
     eta = 0.5
     l1 = 0.0
     l2 = 1E5
     tol = 1E-6

     x_new = 0.0
     lambda = 0.0
     VN = copy(x)

     # Loop while the interval is larger than the allowed tolerance
     while (l2-l1 > tol)

         # Bisection
         lambda = (l2 + l1)/2.0

         # Loop over the elements
         @inbounds for e = 1:ne

             #
             # weigth factor (based on the Optimality criteria)
             #
             #x_new = x[e]
             #if dVol[e]>0.0
             #betae    = max(0.0, -1*dFlex[e])/(lambda*dVol[e])
             betae = -dFlex[e]/(lambda*dVol[e])
             x_new = x[e]*(betae)^eta
             #end

             #
             # Moving limits
             #
             xi = max(1E-3,x[e]*0.9)
             xs = min(1.0, x[e]*1.1)

             VN[e] = max( xi, min( x_new, xs) )

         end #e

         # New volume
         v_new = Volume(ne,VN)

         #
         # New interval
         #
         if (v_new > VL)
             l1 = lambda
         else
             l2 = lambda
         end

     end #while


     # Return the new vector of design variables
     return VN
 end

