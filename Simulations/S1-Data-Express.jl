using Highway
using JLD

Tf = 14400; T = 12; N = 9
t1 = 10800

Lanes = 2
A = fld(Tf-t1, T)+1
Sections = [Highway.Sections1_CR[i] for i in 1:2:length(Highway.Sections1_CR)]
S = length(Sections)

J_in = linspace(0.2/Lanes, 0.6/Lanes, 9); J = length(J_in)

Transition_Diagram = -1*ones(S, Lanes, J)
Speeds_Transition = Array{Float64}(S, Lanes, J)
Densities_Transition = Array{Float64}(S, Lanes, J)

Fluxes_Transition1 = zeros(S, Lanes, J)
Fluxes_Transition2 = zeros(S, Lanes, J)

D_Speeds_Transition = Array{Float64}(S, Lanes, J)
D_Densities_Transition = Array{Float64}(S, Lanes, J)

D_Fluxes_Transition1 = zeros(S, Lanes, J)
D_Fluxes_Transition2 = zeros(S, Lanes, J)

jldopen("S1-ConDatos-Ampli.jld", "w") do file
  for  (num_j, j) in enumerate(J_in)

      v2 = zeros(S, Lanes, A)
      d2 = zeros(S, Lanes, A)
      f2 = zeros(S, Lanes, A)

      f1 = zeros(S, Lanes, A)

      println("Start simulations J_in = $j")
      for n = 1:N

          Sense_1, S1 = Highway.Sense1()

          flujo_local2 = zeros(S, Lanes, A)
          densidad_local2 = zeros(S, Lanes, A)
          velocidad_local_promedio2 = zeros(S, Lanes, A)

          flujo_local1 = zeros(S, Lanes, A)

          println("Allocation of arrays finished, start of simulation $n")

          for t = 0:Tf-1

                Merge_Left_Right!(S1[2], S1[1], 2)
                Merge_Right_Left!(S1[1], S1[2], 1)

                f = j*Lanes/P_ramp_S1[1][1]

                if rand() < j
                  Insert_Vehicle!(S1[1].highway, 2t, f*P_ramp_S1[1][2]/Lanes)
                end
                if rand() < j
                  Insert_Vehicle!(S1[2].highway, 2t+1, f*P_ramp_S1[1][2]/Lanes)
                end

                for (num_ramp, (x0, lramp)) in enumerate(S1_in_ramp_ampli)
                    Ramp!(x0, lramp, f*P_ramp_S1[1+num_ramp][1], f*P_ampli_S1[1+num_ramp][2], S1[1], 1)
                    num_ramp = x0 = lramp = 0
                end

                for k = 1:Lanes
                  AccelerationNoiseS1(S1[k].highway)
                  DecelerationMove(S1[k])
                  if t > t1
                      Measure_Frequencies!(S1[k].highway, div(t-t1+1, T)+1, Sections, T,
                                  slice(flujo_local2, :, k, :), slice(densidad_local2, :, k, :),
                                         slice(flujo_local1, :, k, :),
                                  slice(Transition_Diagram, :, k, :, :), num_j)
                  end
                end
            end

            for t = 1:A, k = 1:Lanes, s = 1:S
                add_simple!(slice(velocidad_local_promedio2, :, k, :), slice(flujo_local2, :, k, :),
                                                        slice(densidad_local2, :, k, :), s, t, 1.)
            end

            for a = 1:A, k = 1:Lanes, s = 1:S

                f2[s, k, a] += flujo_local2[s, k, a]/N
                d2[s, k, a] += densidad_local2[s, k, a]/N
                v2[s, k, a] += velocidad_local_promedio2[s, k, a]/N

                f1[s, k, a] += flujo_local1[s, k, a]/N
            end


            println("ending of simulation $n")
            flujo_local2 = velocidad_local_promedio2 = densidad_local2 = flujo_local1 = 0
        end

        for k = 1:Lanes, s = 1:S

            Fluxes_Transition1[s, k, num_j] = round(mean(Float64[f for f in f1[s, k, 1:end-1]]), 2)
            Fluxes_Transition2[s, k, num_j] = round(mean(Float64[f for f in f2[s, k, 1:end-1]]), 2)

            D_Fluxes_Transition1[s, k, num_j] = round(std(Float64[f for f in f1[s, k, 1:end-1]]), 2)
            D_Fluxes_Transition2[s, k, num_j] = round(std(Float64[f for f in f2[s, k, 1:end-1]]), 2)

            v_promedio = round(mean(Float64[v for v in v2[s, k, 1:end-1]]), 2)
            D_v_promedio = round(std(Float64[v for v in v2[s, k, 1:end-1]]), 2)

            d_promedio = round(mean(Float64[v for v in d2[s, k, 1:end-1]]), 2)
            D_d_promedio = round(std(Float64[v for v in d2[s, k, 1:end-1]]), 2)

            Speeds_Transition[s, k, num_j] = v_promedio
            Densities_Transition[s, k, num_j] = d_promedio

            D_Speeds_Transition[s, k, num_j] = D_v_promedio
            D_Densities_Transition[s, k, num_j] = D_d_promedio

            if Transition_Diagram[s, k, num_j] == -1
              #println((v_promedio)

              if v_promedio >= 4.5
                Transition_Diagram[s, k, num_j] = 2
              elseif v_promedio < 4.5 && v_promedio > 0
                Transition_Diagram[s, k, num_j] = 1
              elseif v_promedio + D_v_promedio >= 4.5
                Transition_Diagram[s, k, num_j] = 3
              end
            end
            #println(Transition_Diagram)
        end

        v_promedio = D_v_promedio = d_promedio = D_d_promedio = 0
        v2 = d2= 0
        println("Ends simulation J_in =  $j")
    end
    println("All simulations over, arrays writing starts")
    write(file, "Fluxes_Transition1", Fluxes_Transition1)
    write(file, "Fluxes_Transition2", Fluxes_Transition2)
    write(file, "D_Fluxes_Transition1", D_Fluxes_Transition1)
    write(file, "D_Fluxes_Transition2", D_Fluxes_Transition2)
    write(file, "Transition_Diagram", Transition_Diagram)
    write(file, "Speeds_Transition", Speeds_Transition)
    write(file, "Densities_Transition", Densities_Transition)
    write(file, "D_Speeds_Transition", D_Speeds_Transition)
    write(file, "D_Densities_Transition", D_Densities_Transition)
end
println("The end")
