function data=update_data(X,Y,q,iter,data_old)
    % Moving window/horizon of q datapoints
    if iter<q
    % UPDATE DATA for estimation
        data.X=[data_old.X, X];
        data.Y=[data_old.Y,Y];
    else
        data.X=[data_old.X(:,2:end), X];
        data.Y=[data_old.Y(:,2:end),Y];
    end
end