clear;
clc;
Vback=0;
is_45=1;
% Zigzag Calculations
n_zig_cell=4; 
zig_width=8;
zig_length=4;
E0 = Vback; % ενέργεια παραμονής ηλεκτρονίων στα άτομα του νανοαγωγού
t = -2.5; % ενέργεια μεταφοράς ηλεκτρονίων στα άτομα του νανοαγωγού στους 2 άξονες
% Δημιουργία του πίνακα Α
A = E0*diag(ones(1,(n_zig_cell*zig_width))) + t*diag(ones(1,(n_zig_cell*zig_width)-1),1)+ t*diag(ones(1,(n_zig_cell*zig_width)-1),-1);
% Δημιουργία του πίνακα Β
b = t * [0 1 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 1 0];
B = kron(eye(zig_width), b); 

%armchair calculations
arm_zig_cell=8; 
arm_width=4; % isxuei pos arm_width_actual=2*arm_width
arm_length=4;
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
%final B matrices
B0=kron(eye(arm_width), B0);
B1=kron(eye(arm_width), B1);
B2=kron(eye(arm_width), B2);
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
    end
    
    Hd=kron(diag(maskA),A_transitional)+kron(diag(maskB0,+1),B0)+kron(diag(maskB0,-1),B0')+kron(diag(maskB1,+1),B1)+kron(diag(maskB1,-1),B1')+kron(diag(maskB2,+1),B2)+kron(diag(maskB2,-1),B2')+kron(diag(maskB3,+1),B3)+kron(diag(maskB3,-1),B3');
    [n, m] = size(H); 
    Hd(1:n, 1:m) = H;
    H=Hd;
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



