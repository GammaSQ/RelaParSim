//Define variables:
Epsilon0=625000/(22468879468420441*%pi);
singcharge=-5//-1.602177*10^(-19);
c=299792458;
maxparticles=10;
timestep=10^-1;
duration=10^-0;
steps=int(duration/timestep)
rmax=1.5*10^-2;
width=10000/c*rmax
maxcharge=7000;
particle_mass=1//1.672621777*10^(-27)

delay=rmax/c*0.1;

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
//    r=vectoR(pars,par);
    //topar=gammsq.*squeeze(1+c^(-2)*(-sum(vels.^2,2).*sum(r.^2,2)+sum((r.*vels),2).^2))
    topar=(gammsq-1).*squeeze(sum(r+vels,2)).^(2)./squeeze(sum(r.^2,2))+squeeze(sum(r.^2,2))
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
    //frompar=scalrpar.^2-squeeze(sum(r.^2,2))*diag(sum(vels.^2,2))))*diag(gammsq)
    
    ////// SOMEWHERE IN HERE BUG!!!!!
    frompar=scalrpar.^(2)./squeeze(sum(r.^2,2)).*(gammsq-1)+squeeze(sum(r.^2,2))
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
function[elecscal]=elecByPar(scalaR,gammsq)
    elecscal=-gammsq./scalaR.*((particle_number-1)*singcharge)
endfunction

//Define field function:
function[elecVec]=elecFieldFunc(q,r,t);
    affected=rmx-c*t;
    fullyAffected=affected+delay*c;
    distance=norm(r(:,[1:3]));
    elecVec=-r;
endfunction

//relativistic correction:
function[ret]=correction(velocities,gamsq)
    velsq=squeeze(sum(velocities.^2,2))
    ret=((1+gamsq.*velsq./c^2.).*sqrt(gamsq))//.^(-1)
endfunction



//simulation function:
function[output]=iteration(timeline, dt)
    adding=ones(particle_number,particle_number)
    simuv=timeline(:,[11:13],dt)*timestep/2+timeline(:,[4:6],dt)
    simu_time=zeros(particle_number,13)// ATTENTION! WRONG NUMBER OF LINES!!!! Just testcase!
    simu_time(:,[1:13])=[squeeze(timeline(:,[1:3],dt)),simuv(:,:),Gammasq(simuv),squeeze(timeline(:,[8:13],dt))]
    //create a True-Matrix with a False-Identity
    pro=zeros(particle_number,particle_number)
    flash=(pro==(pro+eye(particle_number,particle_number)))
    //create the empty matrix to be our current distribution
    current_dist=hypermat([particle_number-1,13,particle_number])
    //Just a unit-matrix, cause we will need it.
    unit=ones(particle_number,particle_number)
    //this matrix will hold every particle-particle-distance.
    distance=zeros(particle_number-1,particle_number-1)
    //this matrix will keep track of which particle already has distances to which particles.
    addoneup=ones(particle_number,1)
    //iterate over all our timesteps, effectively creating a look in the past.
    disp(flash)
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
        flash=flash&(~check);
        //since rows can be filled to a different point and each row can get different number of particles, we need an iteraion:
        for t=1:particle_number
            selection=thistime(t,:);
            //if no particle is to be added in this row, we can move on.
            if or(selection) then
                for i=1:particle_number
                    if selection(i) then
                        //the addoneup(i) tells us, at which place to insert the current selection
                        current_dist(addoneup(i),:,i)=timeline(i,:,dt-flashback)'
                        //and of course the next time we want to overwrite the next value
                        addoneup(i)=(addoneup(i)+1);
                    end;
                end;
            end;
        end;
    end;
    //disp(current_dist)
    //Now some particles will be left, because they have no history in our simulation.
    //We aproximate their history by adding all particles from our exposition.
    if or(flash) then
        for t=1:particle_number
            selection=flash(:,t);
            if or(selection) then
                //this time we add all the remaining particles until one particle's environment is complete
                current_dist([addoneup(t):$],:,t)=timeline(selection,:,1);
            end;
        end;
    end;
    gam=squeeze(current_dist(:,7,:))
    r=vectoR(current_dist(:,[1:3],:),simu_time(:,[1:3]))
    re=vectoRE(r)
    v=vDiff(current_dist(:,[4:6],:),simu_time(:,[4:6]))
    scalmat=scalaRSQ(r,v,gam)
    R=loRezSQ(scalmat,re)
    //A=vecPotByPar(v,scalmat)
    E=elecByPar(scalmat,gam)
    ve=vectoRE(current_dist(:,[4:6],:))
        disp('BLAAAAAAAAARGH!')
    vSq_ov_c=sum(current_dist(:,[4:6],:).^2,2)./c^2
    prefield=v
    lortra=scalaRSQpar(r,simuv,squeeze(current_dist(:,7,:)))
    test=squeeze(1-vSq_ov_c)
    for p= 1:particle_number
        prefield(:,:,p)=diag(lortra(:,p))*(diag(E(:,p).*test(:,p))*squeeze(re(:,:,p))+diag(E(:,p).*vSq_ov_c(:,p))*squeeze(ve(:,:,p)))
    end
    field=squeeze(sum(prefield,1))'

    vecE=zeros(particle_number-1,3,particle_number)
    for i=1:particle_number
        vecE(:,:,i)=diag(E(:,i))*squeeze(re(:,:,i))
    end
    //calculating the correctional value due to velocity of particle (at v+a*t/2)
    cor=diag(correction(simu_time(:,[4:6]),simu_time(:,7)))
    acc = cor*field.*singcharge./particle_mass
    //creating a return matrix
    ret=zeros(particle_number,13)
    accv=acc*timestep
    testacc=simuv+accv
    if or(testacc>c) then
        tempt=timestep
        global(timestep)=timestep/2
        ret=iteration(timeline,dt)
        global(timestep)=tempt
        output=ret
    end
    vel=testacc-accv/2
    ret(:,[1:3])=simu_time(:,[1:3])+vel*timestep
    ret(:,[4:6])=vel
    ret(:,7)=Gammasq(vel)
    ret(:,[8:10])=squeeze(sum(vecE,1))'
    ret(:,[$-2:$])=acc
    output=ret
endfunction

//exposition:
//Creating starting-conditions:
fs=[grand(maxparticles,2,'unf',-rmax,+rmax),grand(maxparticles,1,'unf',0,width),grand(maxparticles,2,'nor',1000,200),grand(maxparticles,1,'nor',40000000,200)];
R=rmax^2;
particles=fs((sum(fs(:,[1:2]).^2,2)<R),:);
pre=size(particles);
particle_number=pre(1);
zer=zeros(particle_number,13);
metapar=hypermat([particle_number-1,7,particle_number]);
particles(:,7)=Gammasq(particles(:,[4:6]))
for i=1:(particle_number);
    //metapar(1,:,i)=particles(i,[1:6]);
    sel=zeros(particles);
    sel(:,:)=particles
    sel(i,:)=[]
    metapar([1:$],:,i)=sel;
end

gam=squeeze(metapar(:,7,:))
r=vectoR(metapar(:,[1:3],:),particles(:,[1:3]))
re=vectoRE(r)
v=vectoR(metapar(:,[4:6],:),particles(:,[4:6]))
// MAJOR BUG!!!!! r sometimes creates/loses particles, no idea why! (v and r work the same, v doesn't have that problem.)
if ~(size(r) == size(v)) then
    disp('BLAHRG!')
end
scalmat=scalaRSQ(r,v,gam)
R=loRezSQ(scalmat,re)
//A=vecPotByPar(v,scalmat)
E=elecByPar(scalmat,gam)
ve=vectoRE(metapar(:,[4:6],:))
vSq_ov_c=sum(metapar(:,[4:6],:).^2,2)./c^2
prefield=v
lortra=scalaRSQpar(r,particles(:,[4:6]),squeeze(metapar(:,7,:)))
test=squeeze(1-vSq_ov_c)
//    prefield(:,:,p)=diag(E(:,p).*test(:,p))*squeeze(re(:,:,p))+diag(E(:,p).*vSq_ov_c(:,p))*squeeze(ve(:,:,p))
for p= 1:particle_number
    prefield(:,:,p)=diag(lortra(:,p))*(diag(E(:,p).*test(:,p))*squeeze(re(:,:,p))+diag(E(:,p).*vSq_ov_c(:,p))*squeeze(ve(:,:,p)))
end
field=squeeze(sum(prefield,1))'

vecE=zeros(particle_number-1,3,particle_number)
for i=1:particle_number
    vecE(:,:,i)=diag(E(:,i))*squeeze(re(:,:,i))
end

cor=diag(correction(particles(:,[4:6]),particles(:,7)))
acc = cor*field.*singcharge./particle_mass//diag(correction(particles(:,[4:6]),particles(:,7)))*field.*singcharge
zer(:,[1:10])=[particles(:,[1:7]),squeeze(sum(vecE,1))']
zer(:,[$-2:$])=acc

saving=hypermat([particle_number, 13, steps])
saving(:,[1:10],1)=[particles(:,[1:3]), particles(:,[4:6])+diag(correction(particles(:,[4:6]),particles(:,7)))*acc.*(timestep/2), zer(:,[7:10])]
saving(:,[11:13],1)=diag(correction(particles(:,[4:6]),particles(:,7)))*acc

//simulation:
for dt=[1:steps] //first step is the exposition!
    disp("ATTENTION!!!!!")
    next=iteration(saving,dt)
    saving(:,:,dt+1)=next
end
disp(saving)