function [ECI2RTN] = ECI2RTN(r_N, v_N)

ur = r_N/norm(r_N);
un = cross(r_N, v_N)/norm(cross(r_N, v_N));
ut = cross(un, ur);

ECI2RTN = [ur'; ut'; un'];
end

