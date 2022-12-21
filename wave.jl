using Plots
using Roots

function u1(x,t,lamNvec,a,L) # Wave functinon in the first region with D-R BCs
    N = length(lamNvec)
    S=zeros(length(x))
    # println(S)
    for n = 2:N
        lamN = lamNvec[n]
        A = sin.(lamN*x)/(L/2- sin.(2*lamN*L)/(4*lamN))
        u0int = pi * sin.(lamN)/(pi^2-lamN^2)
        B = u0int * cos.(lamN*a*t) 
        C = B*A
        S = S + C
    end
    return S
end

function u2(x,t,lam1,lam2,a1,a2,L1,L2,k2)
    N = length(lam2)
    S = 0
    for n = 2:N
        lamN = lam2[n]
        A = lamN*1
        
        B = 0
        for Q = 2:length(lam1)
            # println(Q)
            # println(Q)
            # println(lam1[Q])
            B1 = 4 *lam1[Q]*pi.*sin.(lam1[Q]*L1)/((pi^2 - lam1[Q]^2)*(2*lam1[Q]*L1-sin.(2*lam1[Q]*L1))) 
            # println(B1)
            B2 = a2*lamN*(cos.(a1*lam1[Q]*t)-cos.(a2*lamN*t))
        
            B3 = ((a2*lamN-a1*lam1[Q])*(a1*lam1[Q]+a2*lamN))
            # println(B3)
            
            BS = B1*B2/B3
            # println(B1)
            if isnan(BS)
            else
                B=B+BS
            end
        end
       
        C = sin(lamN*L2)*sin(lamN*a2*t)/k2
        # println(C)
        D = A*B+C
        S = S + D
    end
    return S

end

function B1(L1,t,lamVec1,a1)
    return u1(L1,t,lamVec1,a1,L1)
end


function lamRes1(lam,L,H2)
    return lam * cos(lam*L) + H2 * sin(lam*L)
end

function findLam1(L,H2)
    a = find_zeros(x->lamRes1(x,L,H2),0,50)
    return a
end

function testDR()  
    L1 = 5
    H2 = 1
    a_c = .01
    lamVec = findLam1(L1,H2)


    xvec= 0:1e-3:L1

    tvec = 0:1:1000
    for t in tvec
        u0 = u1(xvec,t,lamVec,a_c,L1)
        # println(u0)
        plot(xvec,u0)
        ylims!((-1,1))
        xlims!((0,L1))
        savefig("DRPulse2.5/test"*string(t)*".png")

    end
end

    
testDR()

exit()

L1 = 5
L2 = 5
H2 = 0.2
k2 = 2
h2=0.2
a1 = .01
lam1 = findLam1(L1,H2)
lam2 = findLam1(L2,H2)
a2 = 0.01   


xvec= 0:1e-3:L2

tvec = 0:10:1000
for t in tvec
    # println(t)
    st=t
    u0=[]
    # u0 = zeros((length(xvec)))
    # println(u0)
    for x in xvec
        push!(u0, u2(x,t,lam1,lam2,a1,a2,L1,L2,k2))
    end
    # println(u0)
    # println(u0)
    plot(xvec,u0)

    # println(length(xvec))
    # println(length(u0))
    # ylims!((-3000,3000))
    xlims!((0,L2))
    # println("Help ")
    savefig("DRPulse3/test"*string(st)*".png")
    # println("Hi")

end



