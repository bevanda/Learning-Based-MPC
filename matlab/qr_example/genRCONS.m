function [outputArg1,outputArg2] = genRCONS(inputArg1,inputArg2)
    %==========================================================================
    % Generate Robust inequality CONStraints
    %==========================================================================

    length_Fw = size(F_w_N, 1);

    Aineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, N*m+m);
    bineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, 1);
    b_crx = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, n);

    L_i = zeros(n, N*m); % width of the state tube R_i 
    KL_i = zeros(m, N*m); % width of the input tube KR_i
    disp('Generating constraints on inputs...');
    d_i = zeros(n,1);
    for ind = 1:N
        %disp(['u ind: ', num2str(ind)]);

        KL_i = K*L_i;
        KL_i(:, (ind-1)*m + (1:m)) = eye(m);

        Aineq((ind-1)*length_Fu + (1:length_Fu),1:N*m) = F_u*KL_i;
        bineq((ind-1)*length_Fu + (1:length_Fu)) = h_u - F_u*K*d_i;
        b_crx((ind-1)*length_Fu + (1:length_Fu),:) = -F_u*K*(A+B*K)^(ind-1);

        L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
        d_i = (A+B*K)*d_i + d_0;
    end

    L_i = zeros(n, N*m);
    disp('Generating constraints on states...');
    d_i = d_0;
    for ind = 1:N
        %disp(['x ind: ', num2str(ind)]);
        L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];

        if ind == 1
            disp('Generating terminal constraints on states...');
            Aineq(N*length_Fu + (1:length_Fw), :) = F_w_N*[L_i zeros(n,m); zeros(m,N*m) eye(m)];
            bineq(N*length_Fu + (1:length_Fw)) = h_w_N - F_w_N*[d_i; zeros(m,1)];
            b_crx(N*length_Fu + (1:length_Fw),:) = -F_w_N*[(A+B*K)^(ind); zeros(m,n)];

        else

            Aineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),1:N*m) = F_x*L_i;
            bineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx)) = h_x - F_x*d_i;
            b_crx(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),:) = -F_x*(A+B*K)^(ind);

        end

        d_i = (A+B*K)*d_i + d_0;
    end

    ind = N;
    L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

