import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "03-tasks",
	async up() {
		await zap.transaction.single/*sql*/`            
            create table 
                tasks
            ( 
                id uuid primary key,
                title text not null,
                done boolean not null,
                "boardId" uuid references boards(id) not null,
                "createdAt" timestamptz not null,
                "updatedAt" timestamptz not null,
				"isArchived" boolean not null
            );
        `;
	},
	async down() {
		await zap.transaction.single/*sql*/`
            drop table tasks;
        `;
	},
};

export default migration;
