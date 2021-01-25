bool Obstacle::remove_circular_constraints(Obstacle *o1, Obstacle *o2) {
    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    auto logger = Utility::create_logger("Obstacle");

    bool mirror = false;
    if (o1->getCoordinates_i2() + 1 == o2->getCoordinates_i1()){
        mirror = true;
        std::swap(o1, o2);
    }

    if (o1->getCoordinates_i1() - 1 == o2->getCoordinates_i2()) {
        bool j_overlap = hasOverlap(o1->getCoordinates_j1(), o1->getCoordinates_j2(), o2->getCoordinates_j1(), o2->getCoordinates_j2());
        bool k_overlap = hasOverlap(o1->getCoordinates_k1(), o1->getCoordinates_k2(), o2->getCoordinates_k1(), o2->getCoordinates_k2());
        if (j_overlap && k_overlap) {
            logger->debug("obstacle {} is facing obstacle {} on the left side. Working on {} left side and on {} right side.", o1->get_name(), o2->get_name(), o1->get_name(), o2->get_name());
// another obstacle (o2) at the left side of o1
            overlap = true;
// calculate coordinates of area which should be removed
// the area is for both obstacle the same only if there are equally long
            size_t o1_x1 = o1->getCoordinates_i1();
            size_t o2_x2 = o2->getCoordinates_i2();

            size_t o1_y1;
            size_t o2_y1;
            Obstacle::calculate_area_index(o1, o2, &o1_y1, &o2_y1, CoordinateAxis::Y, true);

            size_t o1_y2;
            size_t o2_y2;
            Obstacle::calculate_area_index(o1, o2, &o1_y2, &o2_y2, CoordinateAxis::Y, false);

            size_t o1_z1;
            size_t o2_z1;
            Obstacle::calculate_area_index(o1, o2, &o1_z1, &o2_z1, CoordinateAxis::Z, true);

            size_t o1_z2;
            size_t o2_z2;
            Obstacle::calculate_area_index(o1, o2, &o1_z2, &o2_z2, CoordinateAxis::Z, false);

            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}", o1_x1, o1_y1, o1_y2, o1_z1, o1_z2, o1->get_name());
            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}", o2_x2, o2_y1, o2_y2, o2_z1, o2_z2, o2->get_name());

            size_t o1_size_removing_indices = (o1_y2 + 1 - o1_y1) * (o1_z2 + 1 - o1_z1);
            size_t o2_size_removing_indices = (o2_y2 + 1 - o2_y1) * (o2_z2 + 1 - o2_z1);

            std::vector <size_t> o1_new(o1->getSize_obstacleLeft());
            std::vector <size_t> o2_new(o1->getSize_obstacleRight());

            size_t o1_new_size_left = 0;
            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->getObstacleLeft()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->getObstacleLeft()[o1_counter_old];
            }

            size_t o2_new_size_left = 0;
            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x2, o2_y1, o2_z1, Nx, Ny);
            size_t o2_current_index = o2->getObstacleLeft()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->getObstacleLeft()[o2_counter_old];
            }

            size_t o1_current_y = o1_y1;
            size_t o1_current_z = o1_z1;
            size_t o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
            bool o1_end = false;

            size_t o2_current_y = o2_y1;
            size_t o2_current_z = o2_z1;
            size_t o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
            bool o2_end = false;
            for (; o1_counter_old < o1->getSize_obstacleLeft() && o2_counter_old < o2->getSize_obstacleRight() && !o1_end && !o2_end; o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->getObstacleLeft()[o1_counter_old];
                o2_current_index = o2->getObstacleLeft()[o2_counter_old];
                if (o1_current_index != o1_removing_index) {
                    o1_new.push_back(o1_current_index);
                    o1_new_size_left++;
                } else {
                    o1_current_y++;
                    if (o1_current_y > o1_y2) {
                        o1_current_y = o1_y1;
                        o1_current_z++;
                        if (o1_current_z > o1_z2) {
                            o1_end = true;
                        }
                    }
                    o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
                }
                if (o2_current_index != o2_removing_index) {
                    o2_new.push_back(o2_current_index);
                    o2_new_size_left++;
                } else {
                    o2_current_y++;
                    if (o2_current_y > o2_y2) {
                        o2_current_y = o2_y1;
                        o2_current_z++;
                        if (o2_current_z > o2_z2) {
                            o2_end = true;
                        }
                    }
                    o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
                }
            }

            if (!o1_end) {
                for (; o1_counter_old < o1->getSize_obstacleLeft() && o1_current_z < o1_z2; o1_counter_old++) {
                    o1_current_index = o1->getObstacleLeft()[o1_counter_old];
                    if (o1_current_index != o1_removing_index) {
                        o1_new.push_back(o1_current_index);
                        o1_new_size_left++;
                    } else {
                        o1_current_y++;
                        if (o1_current_y > o1_y2) {
                            o1_current_y = o1_y1;
                            o1_current_z++;
                        }
                        o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
                    }
                }
            }

            if (!o2_end) {
                for (; o2_counter_old < o2->getSize_obstacleRight() && o2_current_z < o2_z2; o2_counter_old++) {
                    o2_current_index = o2->getObstacleLeft()[o2_counter_old];
                    if (o2_current_index != o2_removing_index) {
                        o2_new.push_back(o2_current_index);
                        o2_new_size_left++;
                    } else {
                        o2_current_y++;
                        if (o2_current_y > o2_y2) {
                            o2_current_y = o2_y1;
                            o2_current_z++;
                        }
                        o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
                    }
                }
            }

            for (; o1_counter_old < o1->getSize_obstacleLeft(); o1_counter_old++) {
                o1_new.push_back(o1->getObstacleLeft()[o1_counter_old]);
                o1_new_size_left++;
            }
            o1_new.resize(o1_new_size_left);

            logger->debug("new size of obstacle {}: {}", o1->get_name(), o1_new_size_left);
            size_t *o1_new_data = new size_t[o1_new_size_left];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_left, Patch::LEFT);

            for (; o2_counter_old < o2->getSize_obstacleLeft(); o2_counter_old++) {
                o2_new.push_back(o2->getObstacleLeft()[o2_counter_old]);
                o2_new_size_left++;
            }
            o2_new.resize(o2_new_size_left);

            logger->debug("new size of obstacle {}: {}", o2->get_name(), o2_new_size_left);
            size_t *o2_new_data = new size_t[o2_new_size_left];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_left, Patch::LEFT);
        }
    }

    if (mirror){
        std::swap(o1, o2);
    }
}