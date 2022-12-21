using Plots

struct Pulse #Characteristics of a pulse
    amplitude
    wavelength
    position # Position of the LHS of the pulse
    speed
end

function filter(A::Pulse, x1,x2,x,c,t,f)
    if x < x1+c*t # pulse is defined only on interval x1 to x2
        return 0
    elseif x > x2 + c*t
        return 0
    else
        return f(x,t) # Pulse can have any functional definition between the endpoints
    end
end

function waveFunction(A,lambda,c,x1)
    function wave(x,t)
        return A + A * cos(2*pi/lambda*(x - 0.5*lambda -x1 - c*t) ) #centers the pulse within x1 and x2
    end
    return wave
end

function generatePulse(A::Pulse)
    # TODO: Specify location and direction of the pulse
    function pulse(x,t)
        x1 = A.position # Pulse is located left of the origin
        x2 = A.position + A.wavelength       # RHS of the pulse is at the origin
        f = waveFunction(A.amplitude,A.wavelength,A.speed,A.position)
        return filter(A,x1,x2,x,A.speed,t,f)

    end
    return pulse
end

function Reflect(A::Pulse, L,cL)
    c1 = abs(A.speed)     # Characteristic Speed of the Incidient wave
    c2 = abs(cL) # Characteristic speed on the OTHER side of "L"
    ampR = A.amplitude*(c2-c1)/(c2+c1) # Calculate amplitude of the reflected wave (absolute value)
    lamR = A.wavelength # Reflected wave has same wavel. as incident
    posR = L + L -A.position - A.wavelength   #Reflected wave is positioned 
    Rpulse = Pulse(ampR,lamR,posR,-A.speed)
    return Rpulse
end

function Transmit(A::Pulse,L,cL)
    c1 = abs(A.speed) # Velocity of the incident wave
    c2 = abs(cL)      # Speed of target material
    ampT = A.amplitude * 2*c2/(c2+c1)
    lamT = c2 *A.wavelength/c1
    posT = -c2/c1*(-A.position+L+A.wavelength) +L +lamT
    Tpulse = Pulse(ampT,lamT,posT,sign(A.speed)*cL)
    return Tpulse
end

function calculatePulses(inc,c1,c2, L, M)
    Reg1 =[]
    Reg2 =[]
    Reg3 =[]
    # R Right-moving 
    # L Left-moving
    # Number indicates region
    R1 = inc            # Incident Wave
    push!(Reg1,R1)
    L1 = Reflect(inc,L,c2)   # Reflects toward negative infinity
    push!(Reg1,L1)
    R2 = Transmit(inc,L,c2)  # Transmits toward M
    push!(Reg2,R2)
    # L2 = Reflect(R2)    # Internal Reflection
    # R3 = Transmit(R2)   # Transmits toward infinity
    # L1 = Transmit(L2)   # Transmits to negative infinity
    # R2 = Reflect(L2)    # internal Reflection
    # L2 = Reflect(R2)
    # R3 = Transmit(L2)
    # generates:
    # L1, R2
    # L2, R3
    # L1. R2
    # L2, R3

    while abs(R2.amplitude) > 1e-3 #Calculate reflections until amplitude goes below threshold
        L2 = Reflect(R2,M,c1)    # Internal Reflection
        push!(Reg2,L2)
        R3 = Transmit(R2,M,c1)   # Transmits toward infinity
        push!(Reg3,R3)
        L1 = Transmit(L2,L,c1)   # Transmits to negative infinity
        push!(Reg1,L1)
        R2 = Reflect(L2,L,c1)
        push!(Reg2,R2)
    end
    # println(length(Reg1))
    # println(length(Reg2))
    # println(length(Reg3))
    #Transform pulse objects to functions
    for i = 1:length(Reg1)
        Reg1[i] = generatePulse(Reg1[i])
    end
    for i = 1:length(Reg2)
        Reg2[i] = generatePulse(Reg2[i])
    end
    for i = 1:length(Reg3)
        Reg3[i] = generatePulse(Reg3[i])
    end


    return Reg1, Reg2, Reg3
end

function testPulseGeneration()
    c1 = 8.0
    c2 = 2.0
    L  = 1.0
    M  = 3.0
    N  = 5.0
    initialPulse = Pulse(1.0,1,0,c1)
    R1,R2,R3 = calculatePulses(initialPulse,c1,c2, L, M)
    println("Neat")
end

function analyze()


end

function testRegions()
    c1 = 2.0
    c2 = 8.0
    L  = 2.0
    M  = 3.0
    N  = 5.0
    initialPulse = Pulse(0.5,1,0,c1)
    R1,R2,R3 = calculatePulses(initialPulse,c1,c2, L, M)
    additionalPulse=Pulse(-0.5,1,-0.5,c1)
    R1b,R2b,R3b = calculatePulses(additionalPulse,c1,c2, L, M)
    

    # A = generatePulse(initialPulse)
    # push!(R1,A)
    # B = generatePulse(Reflect(initialPulse,L,c2))
    # C = generatePulse(Transmit(initialPulse,L,c2))
    # push!(R1,B)
    # push!(R2,C)
    tvec = 0:0.01:5.0
    x1 = 0.0:0.01:L
    x2 = L:0.01:M
    x3 = M:0.01:N
    for i = 1:length(tvec)
        t = tvec[i]
        # plot(x,A.(x,t),label="I")
        # plot!(x,B.(x ,t),label="R1")
        # plot!(x,C.(x,t),label="T1")
        # u1 = R1[1].(x1,t) + R1[2].(x1,t)
        u1 = sum([R1[i].(x1,t) for i in 1:length(R1)])+sum([R1b[i].(x1,t) for i in 1:length(R1b)])
        u2 = sum([R2[i].(x2,t) for i in 1:length(R2)])+sum([R2b[i].(x2,t) for i in 1:length(R2b)])
        u3 = sum([R3[i].(x3,t) for i in 1:length(R3)])+sum([R3b[i].(x3,t) for i in 1:length(R3b)])
        # plot(x1,A.(x1,t)+B.(x1,t))
        plot(x1,u1)
        plot!(x2,u2)
        plot!(x3,u3)
        # plot!(x2,C.(x2,t))
        # plot!(x3,zeros(length(x3)))
        vline!([L,M])
        ylims!((-2,2))
        if i <10
            savefig("nice4/Thing00"*string(i)*".png")
        elseif i < 100
            savefig("nice4/Thing0"*string(i)*".png")
        else
            savefig("nice4/Thing"*string(i)*".png")
        end

    end


end

function varyPulseSpacing()
    c1 = 2.0
    c2 = 8.0
    L  = 2.0
    M  = 3.0
    N  = 5.0
    initialPulse = Pulse(-0.5,1,0,c1)
    R1,R2,R3 = calculatePulses(initialPulse,c1,c2, L, M)


    xAddl = -1:0.01:1 # Vary Second Pulse location; otherwise pulses are identical
    maxR2 = []
  
    for k = 1:length(xAddl)
        additionalPulse=Pulse(0.5,1,xAddl[k],c1)
        R1b,R2b,R3b = calculatePulses(additionalPulse,c1,c2, L, M)
        tvec = 0:0.01:5.0
        x1 = 0.0:0.01:L
        x2 = L:0.01:M
        x3 = M:0.01:N
        maxtemp=0 # Reset max pulse intensity for the current initial pulse spacing
        for i = 1:length(tvec)
            t = tvec[i]
            # plot(x,A.(x,t),label="I")
            # plot!(x,B.(x ,t),label="R1")
            # plot!(x,C.(x,t),label="T1")
            # u1 = R1[1].(x1,t) + R1[2].(x1,t)
            u1 = sum([R1[i].(x1,t) for i in 1:length(R1)])+sum([R1b[i].(x1,t) for i in 1:length(R1b)])
            u2 = sum([R2[i].(x2,t) for i in 1:length(R2)])+sum([R2b[i].(x2,t) for i in 1:length(R2b)])
            u3 = sum([R3[i].(x3,t) for i in 1:length(R3)])+sum([R3b[i].(x3,t) for i in 1:length(R3b)])
            if maximum(abs.(u2)) > maxtemp
                maxtemp = maximum(abs.(u2))
            end

        end
        push!(maxR2,maxtemp)
    end
    plot(xAddl,maxR2)
    savefig("figs/PulseSpacingReverse.png")
end




function testdalambert()
    amp = 0.5
    wavel=2
    c1 = 2.0
    c2 = 3.0
    L = 1.0 # Position of the first wall, transitioning from c1 to c2
    M = 3.0
    N = 5.0
    initialPulse = Pulse(amp,wavel,-wavel,c1)
    A = generatePulse(initialPulse) # Generate Incident pulse # Amplitude, wavelength, LHS location
    Rpulse = Reflect(initialPulse,L,c2)
    B = generatePulse(Rpulse)
    Tpulse = Transmit(initialPulse,L,c2)
    C = generatePulse(Tpulse)
    D = generatePulse(Reflect(Tpulse,M,c1))
    E = generatePulse(Transmit(Tpulse,M,c1))
    # t = 0 
    tvec = 0.0:0.1:8.0

    # R1 0 to L
    # R2 L to M
    # R3 M to N 

    x = -2.0:0.01:N

    for i in range(1,length(tvec);step=1)
        t = tvec[i]
        plot(x,A.(x,t),label="Incident")
        plot!(x,B.(x ,t),label="Reflected")
        plot!(x,C.(x,t),label="Transmitted")
        # plot!(x,D.(x,t),label="R2")
        # plot!(x,E.(x,t),label="T2")
        vline!([L,M],label=nothing)
        ylims!((-2,2))
        plot!(legend=:bottomleft)
        
        # legend
        savefig("TP/Thing"*string(i)*".png")

    end


end
if abspath(PROGRAM_FILE) ==@__FILE__
    # testdalambert()
    # testRegions()
    # testPulseGeneration()
    varyPulseSpacing()
end