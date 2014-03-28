=begin lip

=pod 

=head1 INTRODUCTION

This is an implementation of the recursive procedure described in section 8.6 of
Computational Geometry in C (2nd Ed.) 
Joseph O'Rourke, Cambridge University Press 1998

Given an Annulus R representing n-1 links of an n-link arm A, with the 
circle C of radius l_n, centered on p. Assuming A can reach the target point p,
we know R∩C is non-empty. The procedure then is:
    (1) ∂R⋂C≠∅?
        Choose one of the two points of intersection t
    (2) R⊇C
        Choose any point t on C, say the point furthest from J_0
where ∂R is the boundary of the annulus R and J_0 is the origin of the circle C

=cut

=pod 

=head2 CODE

=cut 

=pod

=head3 Includes

Note that B<Lip::Pod> is not used directly but, rather, necessary for 
generating the internal documentation. However, it is explicitly use-d
to make sure anyone using this module gets it.

=cut

use v5.10;
use strict;
use warnings;
use Lip::Pod;
use Getopt::Long;
no warnings "recursion";

use POSIX;
use boolean;
use Math::Trig;
use Chart::Gnuplot;
use Data::Dump qw(pp);
use Math::Libm qw(hypot);
use List::Util qw(max sum);

use Circle;
use Annulus;

=pod 

=head3 Constant Declarations

=cut

use constant LINK_LENGTHS => (100,25,50,90);
use constant TARGET_POINT => [74,77];
use constant PI => 4 * atan2(1, 1);
use constant ORIGIN => [0,0];
use constant EPSILON => 10**-7;
use constant OO=>-9;

use constant LINE_TYPE       => 1;
use constant LINE_WIDTH      => 4;
use constant CIRCLE_COLOR    => "black";
use constant LABEL_POINT_TYPE => 26;
use constant LABEL_POINT_SIZE => 0.333;

=pod 

=head3 is_reachable()

Determine if the target point is within the annulus

=cut

sub is_reachable{
    my ($link_lengths,$target)=@_;
    my $annulus=new Annulus({link_lengths=>$link_lengths});
    #Is the target point inside the annulus? If so, return true. Else, return false.
    my $d=(($target->[0] - ORIGIN->[0])**2 + 
           ($target->[1] - ORIGIN->[1])**2);
    if($d - $annulus->inner_radius()**2 > EPSILON && fabs($d - $annulus->inner_radius()**2) > EPSILON){
          if($d - $annulus->outer_radius()**2 < EPSILON && fabs($d - $annulus->outer_radius()**2) > EPSILON){
            return true;    
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}

=pod

=head3 find_configuration()

Given a point p to reach and a list of link lengths specifying the arm, first 
determine if p is reachable. If it is find a configuration using recursion.

=cut

sub find_configuration{
    my ($link_lengths,$target, $iterative)=@_;
    if(is_reachable($link_lengths,$target)){
        my $link_configurations=[];
        my $annulus;
        unless($iterative){
            $annulus=recursive_link_configuration($link_configurations,$link_lengths,$link_lengths,$target);
        }
        else{
            $annulus=iterative_link_configuration($link_lengths,$target);
        }
        return $annulus;
    }
    else{
        die("Not Reachable\n");
    }
}

=pod

=head3 recursive_link_configuration()

=cut

sub recursive_link_configuration{
    my ($link_end_points,$link_lengths_all,$links_reduced,$target)=@_;
    my $number_links=()=@$links_reduced;
    if($number_links==2){#bottom of the recursion
        my $circle_a=new Circle({center=>ORIGIN,radius=>$link_lengths_all->[0]});
        my $circle_b=new Circle({center=>$target,radius=>$link_lengths_all->[1]});
        my $intersections=circle_circle_intersection($circle_a,$circle_b);
        my $t;
        unless($intersections==OO){
            $t=select_intersection($intersections);
        }
        else{
            my $r=deg2rad(int(rand(360)));
            $t=[cos($r)*$circle_a->radius()+$circle_a->center()->[0],sin($r)*$circle_a->radius()+$circle_a->center()->[1]];
        }
        $link_end_points->[0]=$circle_a->center();
        $link_end_points->[1]=$t;
        $link_end_points->[2]=$circle_b->center();
        my $annulus=new Annulus({link_end_points=>$link_end_points,link_lengths=>[$link_lengths_all->[0],$link_lengths_all->[1]]});
        return $annulus;
    }
    else{
        my $links_remaining=()=@$links_reduced;
        my @current_links=@$links_reduced;
        my @links_reduced_current=@$links_reduced[0..$links_remaining-2];
        my $annulus=new Annulus({link_end_points=>$link_end_points,link_lengths=>\@links_reduced_current});
        my $current_circle=new Circle({center=>$target,radius=>$links_reduced->[$links_remaining-1]});
        my $intersections=annulus_circle_intersection($annulus,$current_circle);
        my $target_new;
        unless($intersections==OO){
            $target_new=select_intersection($intersections);
        }
        else{
            my $r=deg2rad(int(rand(360)));
            $target_new=[cos($r)*$current_circle->radius()+$target->[0],sin($r)*$current_circle->radius()+$target->[1]];
        }
        $annulus=recursive_link_configuration($link_end_points,$link_lengths_all,\@links_reduced_current,$target_new);
        push @$link_end_points, $target;
        $annulus=new Annulus({link_end_points=>$link_end_points,link_lengths=>\@current_links});
        return $annulus;        
    }
}

=pod

=head3 iterative_link_configuration()

=cut

sub iterative_link_configuration{
    my ($link_lengths,$target)=@_;
    my @targets=();
    my @link_end_points=();
    my $number_links=()=@$link_lengths;
    my @links_reduced=@$link_lengths[0..$number_links-2];
    my $annulus=new Annulus({link_end_points=>[],link_lengths=>\@links_reduced});
    my $current_circle=new Circle({center=>$target,radius=>$link_lengths->[$number_links-1]});
    my $intersections=annulus_circle_intersection($annulus,$current_circle);
    my $t;
    unless($intersections==OO){
        $t=select_intersection($intersections);
        unshift @targets, $t;
    }
    else{
        my $r=deg2rad(int(rand(360)));
        $t=[cos($r)*$current_circle->radius()+$current_circle->center()->[0],sin($r)*$current_circle->radius()+$current_circle->center()->[1]];
        unshift @targets,$t;
    }  
    push @link_end_points, $target;
    unshift @link_end_points, $t;
    my $number_links_annulus=()=@{$annulus->link_lengths()};
    $number_links--;  
    @links_reduced=@$link_lengths[0..$number_links-2];    
    while ($number_links > 2){
        $annulus=new Annulus({link_end_points=>\@link_end_points,link_lengths=>\@links_reduced});
        $current_circle=new Circle({center=>$targets[0],radius=>$link_lengths->[$number_links-1]});
        $intersections=annulus_circle_intersection($annulus,$current_circle);
        my $t;
        unless($intersections==OO){
            $t=select_intersection($intersections);
            unshift @targets, $t;
        }
        else{
            my $r=deg2rad(int(rand(360)));
            $t=[cos($r)*$current_circle->radius()+$targets[0]->[0],sin($r)*$current_circle->radius()+$targets[0]->[1]];
            unshift @targets,$t;
        }  
        unshift @link_end_points,$t;   
        $number_links_annulus=()=@{$annulus->link_lengths()};   
        $number_links--; 
        @links_reduced=@$link_lengths[0..$number_links-2];  
    }
    if($number_links==2){
        my $circle_a=new Circle({center=>ORIGIN,radius=>$link_lengths->[0]});
        my $circle_b=new Circle({center=>$targets[0],radius=>$link_lengths->[1]});
        my $intersections=circle_circle_intersection($circle_a,$circle_b);
        my $t;
        unless($intersections==OO){
            $t=select_intersection($intersections);
            unshift @targets, $t;
        }
        else{
            my $r=deg2rad(int(rand(360)));
            $t=[cos($r)*$circle_a->radius()+$circle_a->center()->[0],sin($r)*$circle_a->radius()+$circle_a->center()->[1]];
            unshift @targets,$t;
        }
        unshift @link_end_points,$t;
        unshift @link_end_points,$circle_a->center();
    }
    $annulus=new Annulus({link_end_points=>\@link_end_points,link_lengths=>$link_lengths});
    return $annulus;  
}

=pod

=head3 select_intersection()

Given an array of intersection points this subroutine returns a random
intersection point.

=cut

sub select_intersection{
    my($intersections)=@_;
    my $number_intersections=()=@$intersections; 
    my $r=int(rand($number_intersections));
    return $intersections->[$r];
}

=pod

=head3 annulus_circle_intersection()

Returns the points of intersection, if any exist, between the two circles of the
annulus and a third circle.
In two partcular cases we return a single value:
    -In the case that the circle is contained within the annulis we return OO.
    -In the case that there is no intersection we return false.
Otherwise an array reference containing up to four intersection points is returned.

=cut

sub annulus_circle_intersection{
    my($annulus,$circle)=@_;
    my $inner_intersection=circle_circle_intersection($annulus->inner_circle(),$circle);
    my $outer_intersection=circle_circle_intersection($annulus->outer_circle(),$circle);
    if($outer_intersection==OO && $inner_intersection==OO){
        return OO;
    }  
    elsif($outer_intersection==OO && !$inner_intersection){
        return OO;
    }
    elsif(!$inner_intersection && !$outer_intersection){
        die("Circle and annulus are disjoint.\n");
    }
    my @intersections;
    if(ref($inner_intersection) eq "ARRAY"){
        foreach my $t (@$inner_intersection){
            unshift @intersections,$t;
        }
    }
    if(ref($outer_intersection) eq "ARRAY"){
        foreach my $t (@$outer_intersection){  
            unshift @intersections,$t;
        }
    }
    return \@intersections;
}

=pod

=head3 circle_circle_intersection()

Returns the points of intersection, if any exist, between two circles.
In two partcular cases we return a single value:
    -In the case that one circle is contained within another we return -9.
    -In the case that there is no intersection we return false.

=cut

sub circle_circle_intersection{
    my($circle_a,$circle_b)=@_;
    #if either of the circles is simply a point then return false
    #this may happen for the inner circle of the annulus.
    if($circle_a->radius()==0 || $circle_b->radius()==0){
        return false;
    }
    #$dx and $dy are the vertical and horizontal distances between the circle centers.
    my $dx = $circle_b->center()->[0] - $circle_a->center()->[0];
    my $dy = $circle_b->center()->[1] - $circle_a->center()->[1];

    #Determine the straight-line distance between the centers. 
    my $d = hypot($dx,$dy);
    #Check for solvability. 
    if ($d - ($circle_a->radius() + $circle_b->radius()) > EPSILON && fabs($d - ($circle_a->radius() + $circle_b->radius())) > EPSILON){
        #no solution. circles do not intersect. 
        return false;
    }
    if ($d - ($circle_a->radius() - $circle_b->radius()) < EPSILON && fabs($d - ($circle_a->radius() - $circle_b->radius())) > EPSILON){
        #oo solutions. one circle is contained in the other 
        return OO;
    }

    #a is the point where the line through the circle
    #intersection points crosses the line between the circle
    #centers.  

    #Determine the distance from point 0 to a.
    my $a = (($circle_a->radius()*$circle_a->radius()) - ($circle_b->radius()*$circle_b->radius()) + ($d*$d)) / (2.0 * $d) ;

    #Determine the coordinates of a. 
    my $x2 = $circle_a->center()->[0] + ($dx * ($a/$d));
    my $y2 = $circle_a->center()->[1] + ($dy * ($a/$d));

    #Determine the distance from a to either of the
    #intersection points.
    my $h=0;
    my $h_temp=($circle_a->radius()**2) - ($a**2);
    if($h_temp < EPSILON){
        $h=0;
    }
    else{
        $h = sqrt($h_temp);
    }

    #Now determine the offsets of the intersection points from a.
    my $rx = -$dy * ($h/$d);
    my $ry = $dx * ($h/$d);

    #Determine the absolute intersection points. 
    my $xi = $x2 + $rx;
    my $xi_prime = $x2 - $rx;
    my $yi = $y2 + $ry;
    my $yi_prime = $y2 - $ry;
    return [[$xi,$yi],[$xi_prime,$yi_prime]]
}

sub plot_annulus {
    my ($annulus,$target,$index)=@_;
    my $reachability = Chart::Gnuplot->new(
        title   => undef,
        border  => undef,
        legend  => undef,
        xtics   => undef,
        ytics   => undef,
        xlabel  => undef,
        ylabel  => undef,
        tmargin => 0,
        bmargin => 0,
        size    => "square",
        output  => "reachability_$index.eps",
        terminal=> "eps enhanced color font 'Menlo,12' linewidth ".LINE_WIDTH
    );
    $reachability->label(
        pointtype => LABEL_POINT_TYPE,
        pointsize => LABEL_POINT_SIZE,
        text => "O",
        position => "0,0"
    );
    $reachability->label(
        pointtype => LABEL_POINT_TYPE,
        pointsize => LABEL_POINT_SIZE,
        text => "p",
        position => "$target->[0],$target->[1]"
    );
    my $center=$annulus->center();
    my $inner_radius=$annulus->inner_radius();
    my %inner_circle;
    $inner_circle{x} = "$inner_radius*cos(t)+$center->[0]";
    $inner_circle{y} = "$inner_radius*sin(t)+$center->[1]";
    my $inner_circle = Chart::Gnuplot::DataSet->new(
        func     => \%inner_circle,
        color    => CIRCLE_COLOR
    );
    my $outer_radius=$annulus->outer_radius();
    my %outer_circle;
    $outer_circle{x} = "$outer_radius*cos(t)+$center->[0]";
    $outer_circle{y} = "$outer_radius*sin(t)+$center->[1]";
    my $outer_circle = Chart::Gnuplot::DataSet->new(
        func     => \%outer_circle,
        color    => CIRCLE_COLOR,
        linetype => LINE_TYPE
    );
    my @link_end_points=@{$annulus->link_end_points()};
    my $start_point=shift @link_end_points;
    my $end_point=shift @link_end_points;
    while($end_point){
        $reachability->arrow(
                        from     => "$start_point->[0],$start_point->[1]",
                        to       => "$end_point->[0],$end_point->[1]",
                        linetype => LINE_TYPE,
                        head     => "off",
                        color    => CIRCLE_COLOR
                    );
        $reachability->label(
                        pointtype => LABEL_POINT_TYPE,
                        pointsize => LABEL_POINT_SIZE,
                        position  => "$end_point->[0],$end_point->[1]"
        );
        $start_point=$end_point;
        $end_point=shift @link_end_points;
    }
    $reachability->plot2d($inner_circle,$outer_circle);
}

=end lip

=cut

my $iterative;
my $arg_result=GetOptions("iterative"=>\$iterative);

my @links=sort {$b <=> $a} (LINK_LENGTHS);
my $annulus=find_configuration(\@links,TARGET_POINT,$iterative);

plot_annulus($annulus,TARGET_POINT,"all");
say join(" ",map{sprintf("%3.3f",rad2deg($_))} @{$annulus->link_positions});

__END__

=pod

=head1 NAME

reachability.pl

=head1 SYNOPSIS

    perl reachability.pl [--iterative]

=head1 SEE ALSO

This is an implementation of the recursive procedure described on page 329 in
section 8.6 of:
    Computational Geometry in C (2nd Ed.) 
    Joseph O'Rourke, Cambridge University Press 1998

=head1 AUTHOR

Adam Russell, E<lt>arussell@cs.uml.eduE<gt>

=cut
