%% Created by: Andera MASIERO

function I = imenlarger(I0, s)

% s has to be a positive integer



if (size(I0,3) ~= 1) && (size(I0,3) ~= 3)

	error('I0 has to be a m x n x k matrix, where the only admissible values for k are 1 and 3.')

end



I = zeros( s*size(I0,1) , s*size(I0,2), size(I0,3));



for i=1:size(I0,1)

    for j = 1:size(I0,2)

        I( (i-1)*s+1:i*s , (j-1)*s+1:j*s, 1 ) = I0(i,j,1);

		if size(I0,3) == 3

			I( (i-1)*s+1:i*s , (j-1)*s+1:j*s, 2 ) = I0(i,j,2);

			I( (i-1)*s+1:i*s , (j-1)*s+1:j*s, 3 ) = I0(i,j,3);

		end

    end

end

