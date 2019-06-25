function data_out = get_data(X,Y,q,i,data_in)
    % Moving window/horizon of q datapoints
    if i<q
        data = data_in;
        data(:,i+1) = [X; Y];
%         data = [data,[X; Y]];
    else
        data = [data_in(:,2:end),[X; Y]];
    end
    data_out = data;
end