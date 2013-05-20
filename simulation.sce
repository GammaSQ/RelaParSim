//Define variables:
c=299792458;

//Define basic geometrics:
function[vector]=vectoR(pars,par)
    unit=ones(squeeze(pars(:,1,:)))'
    r=pars
    r(:,1,:)=(diag(par(:,1))*unit)'-squeeze(pars(:,1,:))
    r(:,2,:)=(diag(par(:,2))*unit)'-squeeze(pars(:,2,:))
    r(:,3,:)=(diag(par(:,3))*unit)'-squeeze(pars(:,3,:))
    vector=r
endfunction

function[vector]=vectoRE(r)
    dist=(sqrt(sum(r.^2,2)))
    re=r
    re(:,1,:)=[r(:,1,:)./dist]
    re(:,2,:)=[r(:,2,:)./dist]
    re(:,3,:)=[r(:,3,:)./dist]
    vector=re
endfunction

function[velocities]=vDiff(pars,par)
    unit=ones(squeeze(pars(:,1,:)))
    vels=pars
    vels(:,1,:)=unit*diag(par(:,1))-squeeze(pars(:,1,:))
    vels(:,2,:)=unit*diag(par(:,2))-squeeze(pars(:,2,:))
    vels(:,3,:)=unit*diag(par(:,3))-squeeze(pars(:,3,:))
    velocities=vels
endfunction

//scalaR  returns the Lorentz-transformated distances as scalars, no direction!
function[scalmat]=scalaRSQ(r, vels, gammsq)
    noop=and(r==0,2)
    disp(size(r))
    topar(noop)=gammsq(noop)//.*squeeze(sum(r(noop,:)+vels(noop,:),2)).^(2)./squeeze(sum(r(noop,:).^2,2))+squeeze(sum(r(noop,:).^2,2))
    topar(~noop)=0
    scalmat=squeeze(topar)
endfunction

function[scalmat]=scalaRSQpar(r,vels,gammsq)
    //every mul is one coordinate
    mul1=diag(vels(:,1))
    mul2=diag(vels(:,2))
    mul3=diag(vels(:,3))
    //multiplying each particle of all particles with their corresponding mul results in a scalar product of all particles with the right particle
    scalrpar=squeeze(sum(cat(3,squeeze(r(:,1,:))*mul1,squeeze(r(:,2,:))*mul2,squeeze(r(:,3,:))*mul3),3))
    //scalrpar is the scalar product (velocitie*(pos-pos0))

    ////// SOMEWHERE IN HERE BUG!!!!!
    frompar=(scalrpar.^(2)./squeeze(sum(r.^2,2)))*diag((gammsq-1))+squeeze(sum(r.^2,2))
    scalmat=squeeze(frompar)
endfunction

//returns Lorentz-transformed VECTORS!
function[lorentz]=loR(scale,re)
    scalmat=scale.^(1/2)
    ret=zeros(particle_number-1,3,particle_number)
    for rt=1:particle_number
        mul=diag(scalmat(:,rt))
        ret(:,[1:3],rt)=squeeze(re(:,[1:3],rt))*mul
    end
    lorentz=ret
endfunction

//returns Lorentz-transformed SQUARED REZIPROKE vectors (for electrical force)!
function[lorentz]=loRezSQ(scale,re)
    scalmat=scale.^(-1)
    ret=zeros(particle_number-1,3,particle_number)
    for rt=1:particle_number
        mul=diag(scalmat(:,rt))
        ret(:,[1:3],rt)=mul*squeeze(re(:,[1:3],rt))
    end
    lorentz=ret
endfunction

//Define basic Lorentz-function:
function[dillation]=Gammasq(v)
    vsq=squeeze(sum(v.^2,2))
    ret=(1-vsq.*(1/(c^2))).^(-1)//((vsq.*(-1/(c^2)))+1).^(-1)
    dillation=ret;
endfunction

//define vector potential A by particles:
//not implemented, gets defined via electrical field.

//define field by particles:
function[elecscal]=elecByPar(scalaR,charge,gammsq)
    elecscal=-gammsq.*charge./scalaR.*(particle_number-1)
endfunction

//relativistic correction:
function[ret]=correction(velocities,gamsq)
    velsq=squeeze(sum(velocities.^2,2))
    ret=((1+gamsq.*velsq./c^2.).*sqrt(gamsq))//.^(-1)
endfunction

//simulation function:
function[output]=iteration(timeline,timestep,simuv,dt,time,elecFieldFunc)
    nl=argn(2)
    if nl>6 then
        error("Expecting 6 arguments at most")
    end
    if nl<2 then
        error("timeline and timestep-size have to be defined!")
    end
    if nl<3 then
        simuv=timeline(:,[4:6],dt)+timeline(:,[11:13],dt)*(timestep/2)
    end
    if nl<4 then
        dt=1
        time=0
    end
    if nl==4 then
        error("current time and timestep-number have both to be set or unset!")
    end
    if typeof(timeline)~="hypermat" then
        error("Expecting a real-valued hypermat-matrix as timeline, got "+typeof(timeline))
    end
    si=size(timeline)
    if si(3)<dt then
        tmp=timeline
        timeline=zeros(size(1),size(2),dt)
        timeline(:,:,1:si(3))=tmp
    end
    siv=size(simuv)
    if siv(1)~=si(1) then
        error("Need a velocitie for EVERY particle!")
    end
    if si(2)~=13 then
        error("Invalid number of particle-information!")
    end
    if size(timestep,'*')~=1|size(dt,'*')~=1|size(time,'*')~=1 then
        error("Expecting scalar for timestep!")
    end
    adding=ones(particle_number,particle_number)
    //simuv=timeline(:,[4:6],dt)+timeline(:,[11:13],dt)*(timestep/2)
    simu_time=zeros(particle_number,13)// ATTENTION! WRONG NUMBER OF LINES!!!! Just testcase!
    simu_time(:,[1:13])=[squeeze(timeline(:,[1:3],dt)),simuv(:,:),Gammasq(simuv),squeeze(timeline(:,[8:13],dt))]
    //create a True-Matrix with a False-Identity
    pro=zeros(particle_number,particle_number)
    flash=(pro==(pro+eye(particle_number,particle_number)))
    //create the empty matrix to be our current distribution
    current_dist=hypermat([particle_number-1,9,particle_number])
    //Just a unit-matrix, cause we will need it.
    unit=ones(particle_number,particle_number)
    //this matrix will hold every particle-particle-distance.
    distance=zeros(particle_number-1,particle_number-1)
    //this matrix will keep track of which particle already has distances to which particles.
    addoneup=ones(particle_number,1)
    //iterate over all our timesteps, effectively creating a look in the past.
    for flashback=1:dt-1
        //calculating our distances, creating only the upper half of a diagonal matrix! (lower half is transponation of upper half.)
        //rsults in quare matrix with particle_number - 1 cols/rows
        for t=2:particle_number
            predist=ones(particle_number-t+1,3)*diag(timeline(t-1,[1:3],dt-flashback+1))
            distance([t-1:$],t-1)=sum((predist-squeeze(simu_time([t:$],[1:3]))).^2,2).^(1/2)
        end;
        //creating a boolean matrix based on our distances. (we allways chech in the area light could travel during the time past)
        pick=distance<(c*timestep*flashback) & distance~=0//did I really write that shit:((c*timestep*(flashback-1))>=distance & distance>(c*timestep*flashback)) ?
        //create a False-Matrix
        tie=zeros(particle_number,particle_number)
        check=(tie~=tie)
        //identity is false (don't need to calc the particle to itself), upper half defining:
        check([2:$],[1:$-1])=pick
        //lower half is pick', pick is false in lower half, therefore just getting both halfs rsults in an or:
        check([1:$-1],[2:$])=check([1:$-1],[2:$])|pick'
        //check which particles haven't already been maped against eachother:
        thistime=flash&check
        //remove particles of this run from our pool:
        //since rows can be filled to a different point and each row can get different number of particles, we need an iteraion:
        for t=1:particle_number;
            selection=thistime(:,t);
            //if no particle is to be added in this row, we can move on.
            if or(selection) then;
                flash(:,t)=flash(:,t)&(~selection)
                sel=[timeline(selection,[1:7],dt-flashback),parmeta(selection,:)]
                si=size(sel)
                current_dist([addoneup(t):(addoneup(t)+si(1)-1)],:,t)=sel;
                addoneup(t)=addoneup(t)+si(1);
            end;
        end;
    end;
    //Now some particles will be left, because they have no history in our simulation.
    //We aproximate their history by adding all particles from our exposition.
    if or(flash) then
        for t=1:particle_number
            selection=flash(:,t);
            if or(selection) then
                //this time we add all the remaining particles until one particle's environment is complete
                current_dist([addoneup(t):$],:,t)=[squeeze(timeline(selection,[1:7],1)),parmeta(selection,:)]
            end;
        end;
    end;
    gam=squeeze(current_dist(:,7,:))
    charge=squeeze(current_dist(:,8,:))
    r=vectoR(current_dist(:,[1:3],:),simu_time(:,[1:3]))
    re=vectoRE(r)
    v=vDiff(current_dist(:,[4:6],:),simu_time(:,[4:6]))
    scalmat=scalaRSQ(r,v,gam)
    R=loRezSQ(scalmat,re)
    //A=vecPotByPar(v,scalmat)
    E=elecByPar(scalmat,charge,gam)
    ve=vectoRE(current_dist(:,[4:6],:))
    vSq_ov_c=sum(current_dist(:,[4:6],:).^2,2)./c^2
    prefield=v
    lortra=scalaRSQpar(r,simuv,simu_time(:,7))//squeeze(current_dist(:,7,:)))
    test=squeeze(1-vSq_ov_c)
    if nl<6 then
        for p= 1:particle_number
            prefield(:,:,p)=diag(lortra(:,p))*(diag(E(:,p).*test(:,p))*squeeze(re(:,:,p))+diag(E(:,p).*vSq_ov_c(:,p))*squeeze(ve(:,:,p)))
        end
    else
        extra=elecFieldFunc(simu_time,time)
        disp(size(extra))
        for p= 1:particle_number
            prefield(:,:,p)=diag(lortra(:,p))*(diag(E(:,p).*test(:,p))*squeeze(re(:,:,p))+diag(E(:,p).*vSq_ov_c(:,p))*squeeze(ve(:,:,p)))
            prefield(:,:,p)=prefield(:,:,p)+scalaRSQ(extra,simu_time(:,4:6),simu_time(:,7))
        end
    end
    field=squeeze(sum(prefield,1))'

    //calculating the correctional value due to velocity of particle (at v+a*t/2)
    cor=diag(correction(simu_time(:,[4:6]),simu_time(:,7)).*parmeta(:,1)./parmeta(:,2))
    acc = cor*field
    //creating a return matrix
    ret=zeros(particle_number,13)
        vel=simuv+acc*timestep/2
        ret(:,[1:3])=simu_time(:,[1:3])+vel*timestep
        ret(:,[4:6])=vel
        ret(:,7)=Gammasq(vel)
        ret(:,[8:10])=squeeze(sum(vecE,1))'
        ret(:,[($-2):$])=acc
        output=ret
endfunction

function[timeline]=simulation(prime_dist,epsv,timestep,elecFieldFunc)
    nl=argn(2)
    if nl<2 then
        error("Wrong number of minimum arguments!")
    end
    if nl>5 then
        error("Wrong number of argumetns! Can take 5 at most!")
    end
    if nl==2 then
        timestep=epsv/c
    end
    if typeof(prime_dist)=="hypermat" then
        si=size(prime_dist)
        dt=si(3)
    else
        dt=1
    end
    //exposition:
    si=size(prime_dist)
    particle_number=si(1)
    particles=[prime_dist(:,[1:6]),Gammasq(prime_dist(:,[4:6]))]
    parmeta=prime_dist(:,[7:8])
    
    saving=zeros(particle_number, 13, 2)
    saving(:,1:7,1)=particles
    
    time=0
    while time<duration //first step is the exposition!
        disp(time)
        vel=saving(:,[11:13],dt)*timestep/2+saving(:,[4:6],dt)
        if nl==4 then
            next=iteration(saving,timestep,vel,dt,time,elecFieldFunc)
        else
            next=iteration(saving,timestep,vel,dt,time)
        end
        testv=sum((vel).^2,2)
        if or(~isreal(testv))|or(~isreal(sum((saving(:,[11:13],dt)*timestep).^2,2))) then
            prestep=timestep/2
            timestep=prestep(1)
        elseif or(testv<0)|or(testv>c^2)|or(sum((next(:,[11:13])*timestep).^2,2)>epsv^2)
            prestep=min(timestep/2,max(testv/(c^2)),min(sum((next(:,[11:13])*timestep).^2,2).^(1/2)/epsv))
            timestep=prestep(1)
        elseif and(testv<((c^2)/2))&and((sum((next(:,[11:13])*timestep).^2,2))<(epsv^2/4))
            dt=dt+1
            time=time+timestep
            prestep=max(timestep*3,min(min(testv/(c^2)),min(sum((next(:,[11:13])*timestep).^2,2).^(1/2)/epsv)))
            timestep=prestep(1)
        else
            next=iteration(saving,timestep,vel,dt,time)
            dt=dt+1
            time=time+timestep
        end
        saving(:,:,dt)=next
    end
//    timeline=iteration(saving,timestep)
    timeline=saving
endfunction

Epsilon0=625000/(22468879468420441*%pi);
singcharge=0.00001//1.602177*10^(-19);
stepsize=10^(-10);
duration=10^(-9);
steps=int(duration/stepsize)
particle_mass=1.672621777*10^(-27)//Creating starting-conditions:
maxparticles=10;
rmax=1.5*10^-2;
width=10000/c*rmax
Epsilon0=625000/(22468879468420441*%pi);
maxcharge=7000;
delay=rmax/c*0.1;
//,grand(maxparticles,1,'unf',0,width),grand(maxparticles,1,'nor',40000000,200)
fs=[grand(maxparticles,3,'unf',-rmax,+rmax),grand(maxparticles,3,'nor',0,1000)];
R=rmax^2;
prime_dist=fs(sum(fs(:,[1:3]).^2,2)<rmax^2,:)
prime_dist(:,7)=singcharge
prime_dist(:,8)=particle_mass

function[scal]=elecField(charge,r)
    scal=charge*(2*%pi*Epsilon0)*r^(-1)
endfunction

function[elecVec]=elFieldFunc(pars,time)
    r=sum(pars(:,1:3).^2,2)
    si=size(r)
    ret=zeros(si(1),1)
    aff=rmax-(time*c)
    fullaff=delay+rmax-(time*c)
    mat1=r>aff^2
    mat2=r>fullaff^2
    partcharge=cos(delay/(2*%pi)*r.^(1/2)-aff)*maxcharge
    if fullaff>0 then
        ret(mat2,:)=elecField(maxcharge,r(mat2))
    else
        ret(~mat2,:)=elecField(maxcharge,r(~mat2))
    end
    if aff<0 &fullaff>0 then
        ret(r<fullaff^2,:)=elecField(partcharge(r<fullaff^2),r(r<fullaff^2))
    end
    if aff>0 then
        ret(mat1,:)=ret(mat1,:)+elecField(partcharge(mat1),r(mat1))
    else
        ret(~mat1,:)=ret(~mat1,:)-elecField(partcharge(~mat1),r(~mat1))
    end
    elecVec=-diag(ret)*pars(:,1:3)
endfunction

disp(simulation(prime_dist,0.0001,stepsize,elFieldFunc))