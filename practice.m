ranges = zeros(1,301);
i = 200;
ranges(i) = 10;


zero = 151;
%ang = angleBetween(i, zero);
ang = (i - zero) * 0.5;
x = ranges(i) * sind(ang);
y = ranges(i) * cosd(ang);

