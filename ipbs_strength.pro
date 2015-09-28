function ipbs_strength,comp,i,IL,x

xm = x / 1000d0 ;a kG

case IL of
  0: begin ;TRANSITION2
      case comp of
        'N_PI':  case i of
                    0: C = [1.019e-1 ,-4.220e-3,-1.334e-2, 2.952e-3]
                    1: C = [-2.566e-4, 3.714e-2,-1.443e-2, 1.785e-3]
                    2: C = [-1.019e-1, 4.219e-3, 1.334e-2,-2.952e-3]
                    ELSE: return,0
                endcase
        'N_SIG_R':  case i of
                    0: C = [-1.019e-1, 4.219e-3,1.334e-2,-2.952e-3]
                    1: C = [-5.746e-2,-1.458e-2,1.583e-2,-2.904e-3]
                    ELSE: return,0
                endcase
        'N_SIG_B':   case i of
                    0: C = [5.772e-2,-2.256e-2,-1.398e-3,1.118e-3]
                    1: C = [1.019e-1,-4.220e-3,-1.334e-2,2.952e-3]
                    ELSE: return,0
                endcase
        ELSE: return,0
        endcase
  end
  1: begin ;TRANSITION3
      case comp of
        'N_PI':  case i of
                    0: C = [-1.019e-1, 4.220e-3, 1.334e-2,-2.952e-3]
                    1: C = [ 2.606e-4,-3.669e-2, 1.444e-2,-1.789e-3]
                    2: C = [ 1.019e-1,-4.220e-3,-1.334e-2, 2.952e-3]
                    ELSE: return,0
                endcase
        'N_SIG_R':  case i of
                    0: C = [0,0,0,0]
                    1: C = [1.019e-1,-4.220e-3,-1.334e-2, 2.952e-3]
                    2: C = [6.786e-2, 1.435e-2,-1.584e-2, 2.904e-3]
                    ELSE: return,0
                endcase
        'N_SIG_B':   case i of
                    0: C = [-6.812e-2, 2.232e-2, 1.415e-3,-1.117e-3]
                    1: C = [-1.019e-1, 4.220e-3, 1.334e-2,-2.952e-3]
                    2: C = [0,0,0,0]
                    ELSE: return,0
                endcase
      ELSE: return,0
      endcase
  end
  2: begin ;TRANSITION1
      case comp of
      'N_PI':    C = [-1.280   ,-4.635e-4,2.127e-8 ,1.606e-7]
      'N_SIG_R': C = [-1.040e-2,2.326e-4 ,1.577e-5 ,-8.117e-7]
      'N_SIG_B': C = [1.040e-2 ,2.331e-4 ,-1.740e-5,-4.781e-7]
      ELSE: return,0
      endcase
  end
  ELSE: return,0
endcase

return,C[0]*xm +C[1]*xm^2d0 +C[2]*xm^3d0 +C[3]*xm^4d0

end
