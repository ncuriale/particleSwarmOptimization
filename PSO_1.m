clc
clear all
format long g

%constants
r1 = randn();
r2 = randn();
c1 = 1.5;%1.5-2
c2 = 2.0;%2-2.5
dt = 1;
w = 0.4;%0.4-1.44
rho=1;

%initial variables
m=50;%# of particles
n=2;%function dims
x=zeros(m,n);
v=zeros(m,n);
xcon=zeros(m,n,2);%constraint on DVs
vp=v;
xp=x;
for j=1:m
    for i=1:n
        x(j,i)=3*randn();
        v(j,i)=3*randn();
        for k=1:2
            if(k==1); sgn=-1; end;
            if(k==2); sgn=+1; end;
            xcon(j,i,k)=5*sgn;
        end
    end
end
x0=x;
v0=v;

%start arrays
fpi=zeros(m,1);
pi=zeros(m,n);

%iterations
tol=0.001;
res=1;
cnt=0;
while(cnt<100)
    
    %function evaluation
    f=zeros(m,1);
    phi=zeros(m,1);
    for j=1:m
        for i=1:(n-1)
            f(j) = f(j) + 100*(( x(j,i+1) - (x(j,i)^2) )^2) + (1-x(j,i))^2;
        end
    end
    
    %penalty function
    for j=1:m%for all points
        for i=1:n%for all dimensions
            for k=1:2%max/min on design variables
                if(k==1); sgn=-1; end;
                if(k==2); sgn=+1; end;
                phi(j)=phi(j)+(max(0,sgn*(x(j,i)-xcon(j,i,k))))^2;
            end
        end
    end
    f(j)=f(j)+rho*phi(j);%update f 
        
    %update pi
    for j=1:m  
        if(cnt==0)  
            fpi(j)=f(j);
            pi(j,:)=x(j,:);
        else
            if(f(j)<=fpi(j)) 
                fpi(j)=f(j);
                pi(j,:)=x(j,:);
            end
        end
    end
    
    %update pg
    [M,I]=min(f);
    if(cnt==0)
        fpg=f(I);
        pg=x(I,:);
    else
        if(f(I)<=fpg)
            fpg=f(I);
            pg=x(I,:);
        end
    end
    
    %update velocity and position
    for j=1:m
        vp(j,:) = w*v(j,:) + (c1*r1*(pi(j,:)-x(j,:))/dt) + (c2*r2*(pg-x(j,:))/dt);
        xp(j,:) = x(j,:) + vp(j,:)*dt;
        for i=1:n
            for k=1:2
                if(k==1)%low con 
                    if(xp(j,i)<xcon(j,i,k))
                        xp(j,i)=xcon(j,i,k);
                    end;
                elseif(k==2)%high con 
                    if(xp(j,i)>xcon(j,i,k))
                        xp(j,i)=xcon(j,i,k);
                    end;
                end
            end
        end
    end    
    
    x=xp;
    v=vp;
    
    %plot
    figure(1)
    plot([-5,-5],[-5,5],'k')
    hold on
    plot([5,5],[-5,5],'k')
    hold on
    plot([-5,5],[-5,-5],'k')
    hold on
    plot([-5,5],[5,5],'k')
    hold on
    scatter(x(:,1),x(:,2))
    hold on
    grid on
    xlim([-6 6])
    ylim([-6 6])
    drawnow
    pause(1)
    hold off
    cnt=cnt+1;
    rho=1.2*rho;
    %disp(cnt,min(f),max(f))
end

