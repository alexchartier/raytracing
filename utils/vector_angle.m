function a = vector_angle(P1, P2)
a = rad2deg(atan2(norm(cross(P1,P2)),dot(P1,P2))); % Angle in degrees

