function ipbs_shift,comp,i,IL,x

xm = x / 1000d0 ;a kG

case IL of
  0: begin ;TRANSITION2
        case comp of
        'N_PI':  case i of
                    0: C = [8.596,-5.492e-1,-4.028e-1,7.882e-2]
                    1: C = [9.757,-9.627e-1,-4.601e-1,9.984e-2]
                    2: C = [8.517,-5.067e-1,-4.212e-1,8.153e-2]
                    ELSE: return,0
                endcase
        'N_SIG_R':  case i of
                    0: C = [8.477,-4.631e-1,-4.402e-1,8.436e-2]
                    1: C = [9.678,-8.777e-1,-4.968e-1,1.052e-2]
                    ELSE: return,0
                endcase
        'N_SIG_B':   case i of
                    0: C = [9.835,-1.049   ,-4.226e-1,9.430e-2]
                    1: C = [8.596,-5.927e-1,-3.839e-1,7.602e-2]
                    ELSE: return,0
                endcase
        ELSE: return,0
      endcase
  end
  1: begin ;TRANSITION3
      case comp of
        'N_PI':  case i of
                    0: C = [-8.516 ,5.064e-1,4.214e-1,-8.157e-2]
                    1: C = [-11.470,9.645e-1,4.695e-1,-9.970e-2]
                    2: C = [-8.558 ,5.517e-1,4.017e-1,-7.866e-2]
                    ELSE: return,0
                endcase
        'N_SIG_R':  case i of
                    0: C = [-1.603e-2,4.348e-3,6.113e-4,-2.495e-4]
                    1: C = [-8.602   ,5.997e-1,3.807e-1,-7.555e-2]
                    2: C = [-11.550  ,1.059   ,4.188e-1,-9.351e-2]
                    ELSE: return,0
                endcase
        'N_SIG_B':   case i of
                    0: C = [-11.380 ,8.720e-1 ,5.007e-1 ,-1.056e-1]
                    1: C = [-8.474  ,4.600e-1 ,4.415e-1 ,-8.453e-2]
                    2: C = [1.566e-2,-4.191e-1,-5.948e-4,2.399e-3]
                    ELSE: return,0
                endcase
      ELSE: return,0
    endcase
  end
  2: begin ;TRANSITION1
      case comp of
      'N_PI':    C = [1.685,3.073e-2,-1.366e-2,1.567e-3]
      'N_SIG_R': C = [1.685,3.072e-2,-1.366e-2,1.568e-3]
      'N_SIG_B': C = [1.685,3.073e-2,-1.366e-2,1.567e-3]
      ELSE: return,0
    endcase
  end
  ELSE: return,0
endcase

return,C[0]*xm^2d0 +C[1]*xm^3d0 +C[2]*xm^4d0 +C[3]*xm^5d0

end
