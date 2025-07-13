%Author: Georgios-Marios Merkouris 
%Electrical and Computer Engineering Student 
%at Democritus University of Thrace

clear;
clc;
is_45=1;
%apply voltages
Vback=0;
v1=0;
%on how many atoms of the graphene sheet v1 will act
rangeupper=9;
rangelower=1;

% Zigzag part Calculations
n_zig_cell=4; 
zig_width=8;
zig_length=10;
E0 = Vback; % on-site energy of electrons
t = -2.5; % hopping energy of electrons
% Creation of matrix A concerning the zigzag part of the graphene sheet
A = E0*diag(ones(1,(n_zig_cell*zig_width))) + t*diag(ones(1,(n_zig_cell*zig_width)-1),1)+ t*diag(ones(1,(n_zig_cell*zig_width)-1),-1);
% Creation of matrix B concerning the zigzag part of the graphene sheet
b = t * [0 1 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 1 0];
B = kron(eye(zig_width), b); 

% Armchair part calculations
arm_zig_cell=8; 
arm_width=4; % arm_width_actual=2*arm_width due to convention 
arm_length=10;
%A_transitional=E0*diag(ones(1,(n_zig_cell*zig_width))) + t*diag(ones(1,(n_zig_cell*zig_width)-1),1)+ t*diag(ones(1,(n_zig_cell*zig_width )-1),-1);
A_transitional=E0*diag(ones(1,(arm_zig_cell*arm_width))) + t*diag(ones(1,(arm_zig_cell*arm_width)-1),1)+ t*diag(ones(1,(arm_zig_cell*arm_width)-1),-1);
b0=[1 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 1 0 0];
b1=[0 0 0 0;
    0 0 0 0;
    1 0 0 0;
    0 0 0 1];
b3=[1 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 1 0 0];
B0=kron(eye(2),b0);
B1=kron(eye(2),b1);
B2 = zeros(8,8);      
B2(3,2) = 1;          
B2(6,3) = 1;          
B2(7,6) = 1;         
B3=kron(eye(2),b3);
%final B matrices of the armchair part
B0=kron(eye(arm_width), B0);%checked
B1=kron(eye(arm_width), B1);%checked
B2=kron(eye(arm_width), B2);%
% adding additional couplings (e version fix)
n      = size(B2,1);
step   = 8;       % distance between successive rows (and cols)
offset = 3;       % row–col difference (10−7 = 3)
for r = 10:step:n
    c = r - offset;
    if c<=n
        B2(r,c) = 1;
    end
end
B3=kron(eye(arm_width), B3);

% Δημιουργία του Χαμιλτονιανού Πίνακα
if is_45==0
    H = kron(eye(zig_length), A)+ kron(diag(ones(zig_length-1,1), +1), B)+ kron(diag(ones(zig_length-1,1), -1), B');
else
    H = kron(eye(zig_length+1), A)+ kron(diag(ones(zig_length,1), +1), B)+ kron(diag(ones(zig_length,1), -1), B');
    %mask to place matrice A_transitional in H
    maskA = [zeros(zig_length+arm_length,1)];
    maskA(zig_length+2:end) = 1;
    maskB0=[zeros(zig_length+arm_length-1,1)];
    maskB1=[zeros(zig_length+arm_length-1,1)];
    maskB2=[zeros(zig_length+arm_length-1,1)];
    maskB3=[zeros(zig_length+arm_length-1,1)];
    for j=1:arm_length
        if(j<=4)
            if j == 1 %% for b0
                maskB0(zig_length+1)=1;
                B_right=B0;
            elseif mod(j,4) == 2 %% for b1
                if(zig_length+j<=zig_length+arm_length-1)
                maskB1(zig_length+j)=1;
                end
                B_right=B1;
            elseif mod(j,4) == 3 %% for b2
                if(zig_length+j<=zig_length+arm_length-1)
                maskB2(zig_length+j)=1;
                end
                B_right=B2;
            elseif (mod(j,4)==0) && (j>0)%% for b3
                if(zig_length+j<=zig_length+arm_length-1)
                maskB3(zig_length+j)=1;
                end
                B_right=B3;
            end
        % else
        %     if mod(j,4) == 1 %% for b1
        %         if(zig_length+j<=zig_length+arm_length-1)
        %         maskB1(zig_length+j)=1;
        %         end
        %         B_right=B1;
        %     elseif mod(j,4) == 2 %% for b2
        %         if(zig_length+j<=zig_length+arm_length-1)
        %         maskB2(zig_length+j)=1;
        %         end
        %         B_right=B2;
        %     elseif mod(j,4)==3%% for b3
        %         if(zig_length+j<=zig_length+arm_length-1)
        %         maskB3(zig_length+j)=1;
        %         end
        %         B_right=B3;
        %     end        
        % end
            else
        % υπολογίζουμε index από 0->B1, 1->B2, 2->B3
        cycle = mod(j-2, 3);
        pos   = zig_length + j;
        if cycle == 0        %% B1
            if pos <= zig_length+arm_length-1
                maskB1(pos) = 1;
            end
            B_right = B1;
        elseif cycle == 1    %% B2
            if pos <= zig_length+arm_length-1
                maskB2(pos) = 1;
            end
            B_right = B2;
        elseif cycle == 2    %% B3
            if pos <= zig_length+arm_length-1
                maskB3(pos) = 1;
            end
            B_right = B3;
        end
        end
    end
    
    Hd=kron(diag(maskA),A_transitional)+kron(diag(maskB0,+1),B0)+kron(diag(maskB0,-1),B0')+kron(diag(maskB1,+1),B1)+kron(diag(maskB1,-1),B1')+kron(diag(maskB2,+1),B2)+kron(diag(maskB2,-1),B2')+kron(diag(maskB3,+1),B3)+kron(diag(maskB3,-1),B3');
    [n, m] = size(H); 
    Hd(1:n, 1:m) = H;
    H=Hd;
    temp1=zig_width*n_zig_cell*rangelower;
    temp2=rangeupper*zig_width*n_zig_cell;
    for i=temp1+1:temp2
        H(i,i)=H(i,i)+v1;
    end
end
% Δημιουργία των πινάκων Sigma1_template και Sigma2_template
% που είναι μηδενικοί πίνακες εκτος των στοιχείων (1,1) και (4,4)
% που ορίζονται ίσα με 1 και χρησιμοποιούνται για το σχηματισμό των Σ1 και
% Σ2
Sigma1_template=diag([1 zeros(1,zig_length+arm_length-1)]);
Sigma2_template=diag([zeros(1,zig_length+arm_length-1) 1]);

% Δήλωση της μιγαδικής παραμέτρου iη
zplus= 1i*1e-12; 
ii=1;
for EE = -0.5 : 0.01 : +0.5
    E=t*EE;
    % η επαναληπτική μέθοδος προσέγγισης της επιφανειακής συνάρτησης Green
    % για την αριστερή επαφή
    gs1 = inv((E+zplus)*eye((arm_zig_cell*arm_width)) - A); % 1η προσέγγιση της g
    change = 1;
    while change > 1e-6 % συνθήκη που ορίζει την επιθυμητή ακρίβεια σύγκλισης
        Gs = inv((E+zplus)*eye((arm_zig_cell*arm_width)) - A - B' * gs1 * B); % υπολογισμός της προσέγγισης της g 
        % σε κάθε βήμα
        change = sum(sum(abs(Gs - gs1))) / sum(sum(abs(Gs) + abs(gs1))); %υπολογισμός 
        % της απόκλισης της τρέχουσας (Gs) και της προηγούμενης
        % προσέγγισης(gs1)
        gs1 = 0.5 * Gs + 0.5 * gs1; %Ενημερώση της προσέγγισης gs1
        % ως μέσου όρου των δύο τελευταίων για σταθερότερη σύγκλιση
    end
    sig1= B'*gs1*B; % υπολογισμός πίνακα Σ1
    Sigma1=kron(Sigma1_template,sig1); % Σχηματισμός πίνακα Σ1

    % για την δεξιά επαφή
    B=B_right;
    gs2 = inv((E+zplus)*eye((arm_zig_cell*arm_width)) - A); % 1η προσέγγιση της g
    change = 1;
    while change > 1e-6 % συνθήκη που ορίζει την επιθυμητή ακρίβεια σύγκλισης
        Gs = inv((E+zplus)*eye((arm_zig_cell*arm_width)) - A - B * gs2 * B'); % υπολογισμός της προσέγγισης της g 
        % σε κάθε βήμα
        change = sum(sum(abs(Gs - gs2))) / sum(sum(abs(Gs) + abs(gs2))); %υπολογισμός 
        % της απόκλισης της τρέχουσας (Gs) και της προηγούμενης
        % προσέγγισης(gs2)
        gs2 = 0.5 * Gs + 0.5 * gs2; %Ενημερώση της προσέγγισης gs2
        % ως μέσου όρου των δύο τελευταίων για σταθερότερη σύγκλιση
    end
    sig2= B*gs2*B'; % υπολογισμός πίνακα Σ2
    Sigma2=kron(Sigma2_template,sig2); % Σχηματισμός πίνακα Σ2

    Gamma1=1i*(Sigma1-Sigma1'); % Υπολογισμός Πίνακα Γ1
    Gamma2=1i*(Sigma2-Sigma2'); % Υπολογισμός Πίνακα Γ2
    GR=inv(((E+zplus)*eye((n_zig_cell * zig_width * zig_length)+(arm_length*arm_width*arm_zig_cell)))-H-Sigma1-Sigma2); % Υπολογισμός καθυστερημένης συνάρτησης Green
    T(ii)=real(trace(Gamma1*GR*Gamma2*GR')); % Υπολογισμός Κανονικοποιημένης Αγωγιμότητας 
    %για τη συγκεκριμένη ενέργεια Ε
    Eaxis(ii)=EE;
    if EE==0 % για να μη μου κάνει spike η T(E) για EE=0
    T(ii)=T(ii-1);
    end
    ii=ii+1;
end
% Διάγραμμα Ενέργειας - Κανονικοποιημένης Αγωγιμότητας 
figure;
plot(T, Eaxis, 'b-', 'LineWidth', 2.5);
xlabel('Conductance (h/2q^2)G');
ylabel('Energy (E-Ef)/t');
title('Energy vs Conductance Diagram');
grid on;
xlim([0, 10]);

Anum=n_zig_cell*zig_width; % Number of atoms in the column of the Full Hamiltonian
AnumNew=n_zig_cell*arm_length;
Atoms=zeros((Anum*zig_length)+(n_zig_cell*arm_length*arm_width),(Anum*zig_length)+(n_zig_cell*arm_length*arm_width));
Atoms=[ones(Anum*zig_length,Anum*zig_length)];% Sets 1 for each atom in the FULL GNR
%Atoms(1:NWnew*4,zig_length:(zig_length+NLnew))=1;
for i=1:(zig_length+arm_length)*arm_zig_cell*arm_width
    for j=i:(i+(arm_zig_cell*arm_width*arm_length))
        Atoms(i,j)=1;
    end
end
[n, m] = size(Atoms);
for i=1:n
    for j=((zig_length+arm_length)*arm_zig_cell*arm_width)+1:m
        Atoms(i,j)=0;
    end
end


Atoms2=Atoms(1:n,1:n);
figure; 
spy(Atoms2);

