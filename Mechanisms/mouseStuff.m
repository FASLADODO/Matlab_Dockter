function mouseStuff (object, eventdata)
    C = get (gca, 'CurrentPoint');
    fprintf( '(X,Y) = %d,  %d \n', num2str(C(1,1)), num2str(C(1,2)));
end