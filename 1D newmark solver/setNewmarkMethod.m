function setNewmarkMethod(method)
    global betta gamma;
    
    switch method
        case 1
            betta = 1/4;
            gamma = 1/2;
        case 2
            betta = 0;
            gamma = 1/2;
        case 3
            betta = 1/12;
            gamma = 1/2;
        otherwise
            betta = 0.3025;
            gamma = 0.6;
    end
end