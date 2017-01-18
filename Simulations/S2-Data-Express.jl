using Highway
using JLD

Tf = 14400; T = 12; N = 9
t1 = 10800

Lanes = 2; Lanes_t = 3
A = fld(Tf-t1, T)+1
Sections = [Highway.Sections2_CR[i] for i in 1:2:length(Highway.Sections2_CR)]
S = length(Sections)

J_in = linspace(0.2/Lanes, 0.6/Lanes, 9); J = length(J_in)

peeds_Transition = Array{Float64}(S, Lanes, J)
Densities_Transition = Array{Float64}(S, Lanes, J)
Fluxes_Transition = Array{Float64}(S, Lanes, J)

D_Speeds_Transition = Array{Float64}(S, Lanes, J)
D_Densities_Transition = Array{Float64}(S, Lanes, J)
D_Fluxes_Transition = Array{Float64}(S, Lanes, J)

Fluxes_Transition1 = zeros(S, Lanes, J)
Fluxes_Transition2 = zeros(S, Lanes, J)

Speeds_Transition1 = zeros(S, Lanes, J)
Speeds_Transition2 = zeros(S, Lanes, J)

Densities_Transition1 = zeros(S, Lanes, J)
Densities_Transition2 = zeros(S, Lanes, J)

D_Fluxes_Transition1 = zeros(S, Lanes, J)
D_Fluxes_Transition2 = zeros(S, Lanes, J)

D_Speeds_Transition1 = zeros(S, Lanes, J)
D_Speeds_Transition2 = zeros(S, Lanes, J)

D_Densities_Transition1 = zeros(S, Lanes, J)
D_Densities_Transition2 = zeros(S, Lanes, J)

jldopen("S2-Data-Express.jld", "w") do file
  for  (num_j, j) in enumerate(J_in)

    v0 = zeros(S, Lanes, A)
    d0 = zeros(S, Lanes, A)
    f0 = zeros(S, Lanes, A)

    v1 = zeros(S, Lanes, A)
    d1 = zeros(S, Lanes, A)
    f1 = zeros(S, Lanes, A)

    v2 = zeros(S, Lanes, A)
    d2 = zeros(S, Lanes, A)
    f2 = zeros(S, Lanes, A)

    println("Start simulations J_in = $j")
    for n = 1:N

      Sense_2, S2 = Highway.Sense2()

      flujo_local0 = zeros(S, Lanes, A)
      velocidad_local_promedio0 = zeros(S, Lanes, A)

      flujo_local1 = zeros(S, Lanes, A)
      velocidad_local_promedio1 = zeros(S, Lanes, A)

      flujo_local2 = zeros(S, Lanes, A)
      velocidad_local_promedio2 = zeros(S, Lanes, A)

  	  Arr_Speeds0 = Array(Any, S, Lanes, A)
  	  Arr_Speeds1 = Array(Any, S, Lanes, A)
  	  Arr_Speeds2 = Array(Any, S, Lanes, A)

      for t = 1:A, k = 1:Lanes, s = 1:S
        Arr_Speeds0[s, k, t] = Int8[]
        Arr_Speeds1[s, k, t] = Int8[]
        Arr_Speeds2[s, k, t] = Int8[]
      end

      println("Allocation of arrays finished, start of simulation $n")

      for t = 0:Tf-1

        Merge_Left_Right!(S2[2], S2[1], 2)
        Merge_Right_Left!(S2[1], S2[2], 1)

        Merge_Left_Right!(S2[3], S2[2], 3)
        Merge_Right_Left!(S2[2], S2[3], 2)

        f = j*Lanes/P_ramp_S2[1][1]

        if rand() < j
          Insert_Vehicle!(S2[2].highway, 2t+2, f*P_ramp_S2[1][2]/Lanes)
        end
        if rand() < j
          Insert_Vehicle!(S2[3].highway, 2t+1, f*P_ramp_S2[1][2]/Lanes)
        end

        for (num_ramp, (x0, lramp)) in enumerate(S2_in_ramp_ampli)
            Ramp!(x0, lramp, f*P_ramp_S2[1+num_ramp][1], f*P_ramp_S2[1+num_ramp][2], S2[2], 1)
            num_ramp = x0 = lramp = 0
        end

        for k = 1:Lanes_t
          AccelerationNoiseS2(S2[k].highway)
          DecelerationMove(S2[k])
          if t > t1 && k > 1
            Measure_Fluxes!(S1[k].highway, div(t-t1+1, T)+1, Sections, T,
                        slice(flujo_local0, :, k-1, :), slice(flujo_local1, :, k-1, :),
                        slice(flujo_local2, :, k-1, :))
            Measure_Speeds!(S1[k].highway, div(t-t1+1, T)+1, Sections, T,
                        slice(Arr_Speeds0, :, k-1, :), slice(Arr_Speeds1, :, k-1, :),
                        slice(Arr_Speeds2, :, k-1, :))
          end
        end
      end

      ##### Management of Data
      for t = 1:A, k = 1:Lanes, s = 1:S
        velocidad_local_promedio0[s, k, t] = mean(Arr_Speeds0[s, k, t])
        velocidad_local_promedio1[s, k, t] = mean(Arr_Speeds1[s, k, t])
        velocidad_local_promedio2[s, k, t] = mean(Arr_Speeds2[s, k, t])
      end

      Arr_Speeds0 = Arr_Speeds1 = Arr_Speeds2 = 0


	    Measure_Densities!(f0, v0, d0, A, Lanes, S, N, flujo_local0, velocidad_local_promedio0)
	    Measure_Densities!(f1, v1, d1, A, Lanes, S, N, flujo_local1, velocidad_local_promedio1)
	    Measure_Densities!(f2, v2, d2, A, Lanes, S, N, flujo_local2, velocidad_local_promedio2)


      println("ending of simulation $n")
      flujo_local0 = velocidad_local_promedio0 = flujo_local1 = 0
      velocidad_local_promedio1 = flujo_local2 = velocidad_local_promedio2 = 0
    end

  	writing_Arrays!(Lanes, S, num_j, Fluxes_Transition, Fluxes_Transition1, Fluxes_Transition2,
    			 D_Fluxes_Transition, D_Fluxes_Transition1, D_Fluxes_Transition2, f0, f1, f2)

  	writing_Arrays!(Lanes, S, num_j, Speeds_Transition, Speeds_Transition1, Speeds_Transition2,
    			 D_Speeds_Transition, D_Speeds_Transition1, D_Speeds_Transition2, v0, v1, v2)

  	writing_Arrays!(Lanes, S, num_j, Densities_Transition, Densities_Transition1,
    			 Densities_Transition2, D_Densities_Transition, D_Densities_Transition1,
    			 D_Densities_Transition2, d0, d1, d2)

    f0 = f1 = f2 = v0 = v1 = v2 = d0 = d1 = d2= 0
    println("Ends simulation J_in =  $j")
  end
  println("All simulations over, arrays writing starts")
  write(file, "Fluxes_Transition", Fluxes_Transition)
  write(file, "Fluxes_Transition1", Fluxes_Transition1)
  write(file, "Fluxes_Transition2", Fluxes_Transition2)
  write(file, "D_Fluxes_Transition", D_Fluxes_Transition)
  write(file, "D_Fluxes_Transition1", D_Fluxes_Transition1)
  write(file, "D_Fluxes_Transition2", D_Fluxes_Transition2)

  write(file, "Speeds_Transition", Speeds_Transition)
  write(file, "Speeds_Transition1", Speeds_Transition1)
  write(file, "Speeds_Transition2", Speeds_Transition2)
  write(file, "D_Speeds_Transition", D_Speeds_Transition)
  write(file, "D_Speeds_Transition1", D_Speeds_Transition1)
  write(file, "D_Speeds_Transition2", D_Speeds_Transition2)

  write(file, "Densities_Transition", Densities_Transition)
  write(file, "Densities_Transition1", Densities_Transition1)
  write(file, "Densities_Transition2", Densities_Transition2)
  write(file, "D_Densities_Transition", D_Densities_Transition)
  write(file, "D_Densities_Transition1", D_Densities_Transition1)
  write(file, "D_Densities_Transition2", D_Densities_Transition2)
end
println("The end")
