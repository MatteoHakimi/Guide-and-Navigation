close all 
clc, clear

global P_max eta_prop delta_t delta_h Cl Cd S g phi C_s


% Parametri velivolo

b        = 11        ;
S        = 16.2      ;
AR       = b ^ 2 / S ;
P_max    = 172e3     ;
Cd_0     = 0.022     ;
e        = 0.8       ;
eta_prop = 0.82      ;
beta     = 1e-4      ;
C_s      = 0         ;
g        = 9.81      ;


% Calcolo delle condizioni di trim

W_in             = 14100                                      ;
h_0              = 1000                                       ;
phi              = 0                                          ;
Cl               = sqrt ( Cd_0 * pi * AR * e )                ;
Cd               = Cd_0 + Cl ^2 /( pi *AR * e )               ;
[~, ~, ~, rho_0] = atmosisa( h_0 )                            ;
V_1              = sqrt( 2 * W_in / ( rho_0 * S * Cl ) )      ;
D_0              = 1/2 * rho_0 * S * Cd * V_1 ^2              ;
delta_h          = exp( -beta * h_0 )                         ;
delta_t          = D_0 * V_1 / ( eta_prop * P_max * delta_h ) ;


% Condizioni iniziali

V_0      = V_1       ;
psi_0    = 0         ;
gamma_0  = 0         ;
x_0      = 0         ;
y_0      = 0         ;
z_0      = h_0 + 100 ;          % perturbazione su quota iniziale
W_0      = W_in      ;


% Integrazione delle equazioni con Runge Kutta 4° ordine
dt               = 0.001                                             ;
t_span           = 0 : dt : 100                                      ;
X_0              = [ V_0 , psi_0 , gamma_0 , x_0 , y_0 , z_0 , W_0 ] ;
[t,X]            = ode45( 'Cessna_182t_v1', t_span , X_0 )              ;


%Confronto periodo del fugoide con modello di Lanchester

n = input( 'Si vuole calcolare il periodo del Fugoide? [Y/N]:' , 's' ) ;


if isempty(n)
    n = 'Y' ;
end

if n == 'Y'
    n = 1 ;
else
    n = 0 ;
end


switch n
    
    case 1
        [ ~ , locs ]    = findpeaks ( X ( : , 1 ) )    ;
        T_comp          = mean ( diff ( locs )  * dt ) ;
        omega_Lanch = sqrt ( 2 ) *  g / V_0            ;
        T_Lanch     = 2 * pi / omega_Lanch             ;
        err_perc    = abs ( T_comp - T_Lanch ) * 100   ;
        
    case 0
        
end

%% Plot

% Velocità

figure ( 1 )
plot ( t , X ( : , 1 ) , 'b' )
title ( ' Andamento velocità funzione del tempo t ' )
xlabel ( ' t [ s ] ' )                 
ylabel ( ' V [ m/s ] ' )

% Angolo di rotta

figure ( 2 )
plot ( t , X ( : , 2 ) * 180 / pi , 'b' )
title ( ' Andamento angolo di rotta \psi funzione del tempo t ' )
xlabel ( ' t [ s ] ' )
ylabel ( ' \psi [ Deg ] ')

% Angolo di rampa

figure ( 3 )
plot ( t , X ( : , 3 ) * 180 / pi , 'b' )
title ( ' Andamento angolo di rampa \gamma funzione del tempo t ' )
xlabel ( ' t [ s ] ' )
ylabel ( ' \gamma [ Deg ] ' )

% Quota

figure ( 4 )
plot ( t , X ( : , 6 ) , 'b' )
title ( ' Andamento quota h funzione del tempo t ' ) 
xlabel ( ' t [ s ] ' )
ylabel ( ' h [ m ] ' )

% Massa velivolo

figure ( 5 )
plot( t , X ( : , 7 ) / g , 'b')
title( ' Andamento massa m del velivolo funzione del tempo t ' )
xlabel( ' t [ s ] ' )
ylabel ( ' m [ kg ] ' )

% Traiettoria piano fisico xz

figure ( 6 )
plot ( X( : , 4 ) , X( : , 6 ) , 'b' )
title ( ' Traiettoria nel piano xz ' )
xlabel( ' x [ m ] ' )
ylabel ( ' h [ m ] ' )