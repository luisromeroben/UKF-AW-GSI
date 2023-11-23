function err = rmse(v1,v2)
    err = 0;
    for i=1:length(v1)
        err = err + (v1(i)-v2(i))^2;
    end
    err = sqrt((1/length(v1))*err);
end