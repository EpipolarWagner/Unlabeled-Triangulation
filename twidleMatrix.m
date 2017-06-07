function twidleA = twidleMatrix( a )
twidleA=[
  a(1,1)^2,  a(1,1)*a(2,1),  a(1,1)*a(3,1),  a(2,1)^2,  a(2,1)*a(3,1),  a(3,1)^2
  2*a(1,1)*a(1,2),  a(1,2)*a(2,1)+a(1,1)*a(2,2),  a(1,2)*a(3,1)+a(1,1)*a(3,2),  2*a(2,1)*a(2,2),  a(2,2)*a(3,1)+a(2,1)*a(3,2),  2*a(3,1)*a(3,2)
  2*a(1,1)*a(1,3),  a(1,3)*a(2,1)+a(1,1)*a(2,3),  a(1,3)*a(3,1)+a(1,1)*a(3,3),  2*a(2,1)*a(2,3),  a(2,3)*a(3,1)+a(2,1)*a(3,3),  2*a(3,1)*a(3,3)
  2*a(1,1)*a(1,4),  a(1,4)*a(2,1)+a(1,1)*a(2,4),  a(1,4)*a(3,1)+a(1,1)*a(3,4),  2*a(2,1)*a(2,4),  a(2,4)*a(3,1)+a(2,1)*a(3,4),  2*a(3,1)*a(3,4)
  a(1,2)^2,  a(1,2)*a(2,2),  a(1,2)*a(3,2),  a(2,2)^2,  a(2,2)*a(3,2),  a(3,2)^2
  2*a(1,2)*a(1,3),  a(1,3)*a(2,2)+a(1,2)*a(2,3),  a(1,3)*a(3,2)+a(1,2)*a(3,3),  2*a(2,2)*a(2,3),  a(2,3)*a(3,2)+a(2,2)*a(3,3),  2*a(3,2)*a(3,3)
  2*a(1,2)*a(1,4),  a(1,4)*a(2,2)+a(1,2)*a(2,4),  a(1,4)*a(3,2)+a(1,2)*a(3,4),  2*a(2,2)*a(2,4),  a(2,4)*a(3,2)+a(2,2)*a(3,4),  2*a(3,2)*a(3,4)
  a(1,3)^2,  a(1,3)*a(2,3),  a(1,3)*a(3,3),  a(2,3)^2,  a(2,3)*a(3,3),  a(3,3)^2
  2*a(1,3)*a(1,4),  a(1,4)*a(2,3)+a(1,3)*a(2,4),  a(1,4)*a(3,3)+a(1,3)*a(3,4),  2*a(2,3)*a(2,4),  a(2,4)*a(3,3)+a(2,3)*a(3,4),  2*a(3,3)*a(3,4)
  a(1,4)^2,  a(1,4)*a(2,4),  a(1,4)*a(3,4),  a(2,4)^2,  a(2,4)*a(3,4),  a(3,4)^2];
end
