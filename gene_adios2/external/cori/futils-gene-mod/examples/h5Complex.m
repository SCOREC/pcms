function z = h5Complex(file, dset)
    data = hdf5read(file, dset);
    dims = size(data);
    rank = size(dims,2);
    switch rank
      case {1}
        for i=1:dims(1)
              z(i)=complex(cell2mat(data(i,1).Data(1)), cell2mat(data(i,1).Data(2)));
        end
      case {2}
        for i=1:dims(1)
            for j=1:dims(2)
                z(i,j)=complex(cell2mat(data(i,j).Data(1)), cell2mat(data(i,j).Data(2)));
            end
        end
      case {3}
        for i=1:dims(1)
            for j=1:dims(2)
                for k=1:dims(3)
                    z(i,j,k)=complex(cell2mat(data(i,j,k).Data(1)), cell2mat(data(i,j,k).Data(2)));
                end
            end

        end
    end
