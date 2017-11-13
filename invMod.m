function [Inv] = invMod (num, modulo)
    [G, C, ~] = gcd(num, modulo);
    if G == 1  % The inverse of a(mod b) exists only if gcd(a,b)=1
        Inv = mod(C, modulo);
    else
        disp('Modular multiplicative inverse does not exist for these values')
    end
end